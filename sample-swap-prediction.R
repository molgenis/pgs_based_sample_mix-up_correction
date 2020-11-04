#!/usr/bin/env Rscript

#############################################################
## c.a.warmerdam@umcg.nl - April 2020
## Processes plink polygenic scores and phenotypes and
## compares these between samples to quantify the
## to what degree samples are swapped with each other.
#############################################################

##############################
# Load libraries
##############################
library(MASS)
library(tidyverse)
library(argparse)
library(data.table)
library(ROCR)

##############################
# Define argument parser
##############################

parser <- ArgumentParser(description='')
parser$add_argument('--trait-gwas-mapping',
                    help='path to a tab-delimited file that maps traits to the GWAS summary statistics')
parser$add_argument('--base-pgs-path',
                    help="path to a directory containing polygenic scores. Folder structure should be '<base-pgs-path>/<name-of-gwas-summary-statistic>/full.UGLI.pgs.profile'")
parser$add_argument('--phenotypes-file',
                    help='path to a tab-delimited file holding all processed phenotype data.')
parser$add_argument('--sample-coupling-file', required = FALSE,
                    help=paste0('file containing genotype sample ids in the first column',
                    'and phenotype sample ids in the second column'))
parser$add_argument('--llr-bayes-method', required = TRUE, nargs = "+",
                    help = paste0('the naive bayes method to use for calculating log likelihood ratios.',
                                  '<ewi-discretization | efi-discretization | gaussian> [(average) number of samples per bin]'))
parser$add_argument('--out',
                    help='path to output directory')

##############################
# Define functions
##############################

plotSigmoid <- function(estimate, actual, covariates, logitModel) {
  
  covariateLabels <- colnames(covariates)
  
  newData <- data.frame(
    estimate=seq(floor(min(estimate)), ceiling(max(estimate)), 0.05))

  for (lab in covariateLabels) {

    if (is.factor(covariates %>% pull(lab))) {
      newData[, lab] <- sample_n(covariates[, lab], nrow(newData), replace = T)
    }
    else {
      newData[, lab] <- mean(covariates %>% pull(lab))
    }
  }
  
  dat <- data.frame(estimate, actual)
  
  h <-  data.frame(estimate, actual) %>% group_by(actual) %>%
    mutate(breaks = cut(estimate, 
                        breaks=seq(floor(min(estimate)) - 0.025, ceiling(max(estimate)) + 0.025, 0.05), 
                        labels=seq(floor(min(estimate)), ceiling(max(estimate)), 0.05), 
                        include.lowest=TRUE),
           breaks = as.numeric(as.character(breaks))) %>%
    group_by(actual, breaks) %>% 
    summarise(n = n()) %>%
    mutate(pct = ifelse(actual==0, n/sum(n), 1 - n/sum(n))) 
  
  sigmoid <- ggplot() +
    geom_segment(data=h, size=4, show.legend=FALSE,
                 aes(x=breaks, xend=breaks, y=actual, yend=pct, colour=factor(actual))) +
    geom_segment(dat=dat[dat$actual==0,], aes(x=estimate, xend=estimate, y=0, yend=-0.04), size=0.2, colour="grey30") +
    geom_segment(dat=dat[dat$actual==1,], aes(x=estimate, xend=estimate, y=1, yend=1.04), size=0.2, colour="grey30") +
    geom_line(data=data.frame(x=seq(floor(min(estimate)), ceiling(max(estimate)), 0.05), 
                              y=predict(logitModel, 
                                        newdata=newData,
                                        type="response")), 
              aes(x,y), colour="grey50", lwd=1) +
    theme_bw(base_size=12)
  print(sigmoid)
}

# Function that takes (a set of) estimates, actual values, confounders, 
# and a linear model.

# The function predicts actual values from the independent variables 
# (estimates and confounders) using the given model, 
# and returns the deviance of these predicted values
# to the actual values.
devianceFromOlsRegressionLine <- function(estimate, actual, covariates, olsModel) {
  # Predict actual values using the given model and the supplied independent variables
  predictedActual <- predict(olsModel,
                             cbind(estimate, covariates), 
                             type = "response")
  
  # Calculate the residual
  deviance <- actual - predictedActual

  return(unname(deviance))
}

# Function that takes (a set of) estimates, actual values, covariates, 
# and a logistic regression model.

# The function predicts actual values from the independent variables 
# (estimates and covariates) using the given model, 
# and returns the deviance of these predicted values
# to the actual values.
devianceFromLogitRegressionLine <- function(estimate, actual, covariates, logitModel) {
  # Predict actual values using the given model and the supplied independent variables.
  predictedActual <- predict(logitModel, 
                             cbind(estimate, covariates), 
                             type = "response")
  
  # Calculate the distance between the actual values and the predicted actual values.
  distanceBetweenPredictedAndActual <- actual - predictedActual
  
  # Calculate the deviance residuals from the distance between the actual values and
  # the predicted actual values.
  deviance <- sqrt(2 * -log(1 - abs(distanceBetweenPredictedAndActual)))
  deviance <- sign(distanceBetweenPredictedAndActual) * deviance
  
  return(unname(deviance))
}

# Function that takes (a set of) estimates, actual values, covariates, 
# and an ordered logistic regression model.

# The function predicts actual values from the independent variables 
# (estimates and covariates) using the given model, 
# and returns the deviance of these predicted values
# to the actual values.
devianceFromOrderedLogitModel <- function(estimate, actual, covariates, orderedLogitModel) {
  # Predict actual values using the given model and the supplied independent variables.
  predictedProbabilities <- predict(orderedLogitModel, 
                                    cbind(estimate, covariates), 
                                    type = "probs")

  # Determine the predicted probability for the actual class
  predictedProbabilityForActualClass <- predictedProbabilities[
    cbind(seq_len(length(actual)), actual)]
  
  # Set the probabilities for the actual class to 0
  predictedProbabilities[
    cbind(seq_len(length(actual)), actual)] <- 0
  
  # Get the column for which the probability is largest,
  # excluding the column corresponding to the actual class
  competingProbability <- max.col(predictedProbabilities, "first")
  # Get the corresponding classes as an ordered factor
  competingClass <- factor(levels(actual)[competingProbability], 
                           levels = levels(actual), ordered = T)

  # Calculate the deviance residuals from the probability for the actual class
  deviance <- sqrt(2 * -log(predictedProbabilityForActualClass))
  # Add direction:
  # When the probabilities gravitate towards > actual class: negative direction
  # When the probabilities gravitate towards < actual class: positive direction
  deviance <- ((actual > competingClass) - (competingClass > actual)) * deviance
  
  return(deviance)
}

# Function that returns a function for calculating residuals.
# The model determining residuals is dependent on the type of response data type.
residualsFunConstructor <- function(estimate, actual, covariates, responseDataType = "continuous") {
  
  # Return a linear model in case 'actual', is a continuous data type.
  if (responseDataType == "continuous") {
    olsModel <- lm(actual ~ estimate + . + .^2, data = covariates)
    
    print("R-squared:")
    print(summary(olsModel))
    
    # Define residuals function based on the linear regression model.
    residualsFun <- function(estimate, actual, covariates) {
      return(devianceFromOlsRegressionLine(estimate = estimate, 
                                           actual = actual, 
                                           covariates = covariates, 
                                           olsModel = olsModel))
    }
    
    return(residualsFun)
    
  # Return a logistic model in case 'actual', is a binary data type.
  } else if (responseDataType == "binary") {
    logitModel <- glm(actual ~ estimate + . + .^2, 
                      family = binomial(link='logit'),
                      data = covariates)
    
    print(summary(logitModel))
    
    plotSigmoid(estimate = estimate, actual = actual, covariates = covariates, logitModel = logitModel)
    
    # Define residuals function using the logistic regression model.
    residualsFun <- function(estimate, actual, covariates) {
      return(devianceFromLogitRegressionLine(estimate = estimate, 
                                             actual = actual, 
                                             covariates = covariates, 
                                             logitModel = logitModel))
    }
    
    return(residualsFun)
    
  # Return residuals function based on ordered logit in case 'actual' is an ordinal data type,
  # other than binary
  } else if (responseDataType == "ordinal") {
    classes <- sort(unique(actual))
    actualOrdered <- factor(actual, levels = classes, ordered = TRUE)
    
    orderedLogitModel <- polr(actualOrdered ~ estimate + . + .^2, 
                              method = "logistic", 
                              data = covariates)
    
    print(summary(orderedLogitModel))
    
    # Define residuals function using the ordered logistic regression model.
    residualsFun <- function(estimate, actual, covariates) {
      actualOrdered <- factor(actual, levels = classes, ordered = TRUE)
      
      return(devianceFromOrderedLogitModel(estimate = estimate,
                                           actual = actualOrdered,
                                           covariates = covariates,
                                           orderedLogitModel = orderedLogitModel))
    }
    
    return(residualsFun)
  }
}


# Define function for calculating scaled residuals.
# Output is a matrix with rows representing actual values and columns representing expected values
calculate.scaledResiduals <- function(estimate, actual, covariates, responseDataType = "continuous") {
  # Get a function that returns the deviance of an observation to the regression line.
  # This will either be the 'deviance residuals' from a logistic model when the response data type is binary,
  # or this will be the regular residuals corresponding to a linear model.
  residualsFun <- residualsFunConstructor(estimate = estimate,
                                          actual = actual,
                                          covariates = covariates,
                                          responseDataType = responseDataType)

  # extract residuals.
  residuals <- residualsFun(estimate = estimate, 
                            actual = actual, 
                            covariates = covariates)

  # Get the mean and standard deviation to express residuals as z-score-like values.
  residuals.mean <- mean(residuals)
  residuals.sd <- sd(residuals)
  rm(residuals)

  residualsMatrix <- sapply(estimate, function(estimateValue) {

    return(residualsFun(estimate = estimateValue, 
                 actual = actual, 
                 covariates = covariates))
  })
  
  #residualsMatrix <- (residualsMatrix - residuals.mean) / residuals.sd
  return(residualsMatrix)
}

adaptedEqualWidthIntervals <- function(x, nBins, minFrequencyInTails) {
  # Get breaks for outermost null residuals, so that the outermost N null residuals
  # on each tail are represented in their respective bin.
  n <- length(x)
  
  if (n >= (minFrequencyInTails * 2)) {
    
    upperTailLowerBound <- sort(x, partial = n - minFrequencyInTails)[n - minFrequencyInTails]
    upperTailUpperBound <- max(x) + 1
    
    lowerTailUpperBound <- sort(x, partial = minFrequencyInTails)[minFrequencyInTails]
    lowerTailLowerBound <- min(x) - 1
    
    # Separate the space within the outermost bins into a nBins - 2
    binSize <- (upperTailLowerBound - lowerTailUpperBound) / (nBins - 2)
    breaks <- c(lowerTailLowerBound, lowerTailUpperBound,
                lowerTailUpperBound + binSize * (1:(nBins - 3)), 
                upperTailLowerBound, upperTailUpperBound)
  } else {
    
    breaks <- c(min(x) - 1, max(x) + 1)
  }
  
  return(unique(breaks))
}

# Function for converting a matrix of scaled residuals to log likelihood ratios.
# Bins are selected while maintaining a regular bandwidth
scaledResidualsFilteredToLlr.naiveBayes.evenWidthBins <- function(
  scaledResiduals, samplesToFitInGaussian = T, 
  averageSamplesPerBin = 25, minFrequencyInTails = 10) {
  
  # Extract the null-residuals;
  # the residuals belonging to the matches that are assumed to be correct.
  nullResiduals <- scaledResiduals[lower.tri(scaledResiduals, diag = TRUE) 
                                   & upper.tri(scaledResiduals, diag = TRUE)
                                   & samplesToFitInGaussian]
  
  # Extract the alternative residuals; 
  # the residuals belonging to the matches that are assumed to be sample-swaps.
  alternativeResiduals <- scaledResiduals[(lower.tri(scaledResiduals) | upper.tri(scaledResiduals)) 
                                          & samplesToFitInGaussian]

  scaledResidualsFiltered <- scaledResiduals[samplesToFitInGaussian, ]
  
  message("Extracted scaled residuals")
  
  rm(scaledResiduals)
  gc()
  
  # Store the dimensions of the scaled residual matrix, as well as rownames and colnames.
  returnDimensions <- dim(scaledResidualsFiltered)
  returnColumnNames <- colnames(scaledResidualsFiltered)
  returnRowNames <- rownames(scaledResidualsFiltered)
  
  # Calculate the number of bins given the number of null samples per bin and
  # the total number of null residuals
  nBins <- as.integer(length(nullResiduals) / averageSamplesPerBin)

  message("Obtaining breaks")
  
  breaks <- adaptedEqualWidthIntervals(nullResiduals, nBins, minFrequencyInTails)
  breaks[1] <- min(scaledResidualsFiltered) - 1
  breaks[length(breaks)] <- max(scaledResidualsFiltered) + 1
  
  message("Breaks created")

  nullTiles <- cut(nullResiduals, breaks = breaks, labels = FALSE)
  alternativeTiles <- cut(alternativeResiduals, breaks = breaks, labels = FALSE)
  
  message("Breaks applied")

  # Get the density / likelihood of the null residuals for each of the bins.
  nullLikelihoods <- sapply(
    1:nBins,
    function(bin) sum(nullTiles == bin) / length(nullTiles))
  
  # Get the density / likelihood of the alternative residuals for each of the bins.
  alternativeLikelihoods <- sapply(
    1:nBins,
    function(bin) sum(alternativeTiles == bin) / length(alternativeTiles))
  
  message("Alternative likelihoods calculated")
  
  message("Calculated likelihoods")

  # Remove the null and alternative tiles to clear memory.
  rm(nullTiles)
  rm(alternativeTiles)
  gc()
  
  likelihoodRatios <- sapply(scaledResidualsFiltered, function(residual) {
    return(alternativeLikelihoods[cut(residual, breaks = breaks, labels = F)] / 
             nullLikelihoods[cut(residual, breaks = breaks, labels = F)])
  })
  
  message("Calculated likelihood ratios")

  # Apply the same dimensions, row names and column names as where used for the input matrix.
  dim(likelihoodRatios) <- returnDimensions
  rownames(likelihoodRatios) <- returnRowNames
  colnames(likelihoodRatios) <- returnColumnNames
  
  # Return the log-transformed likelihood ratios.
  return(log(likelihoodRatios))
}

scaledResidualsToLlr.naiveBayes.evenWidthBins <- function(
  scaledResiduals, actual, responseDataType = "continuous", 
  averageSamplesPerBin = 25, minFrequencyInTails = 10) {
  
  # Store the dimensions of the scaled residual matrix, as well as rownames and colnames.
  logLikelihoodRatios <- matrix(nrow = nrow(scaledResiduals), 
                                ncol = ncol(scaledResiduals),
                                dimnames = list(rownames(scaledResiduals), colnames(scaledResiduals)))
  
  # Check whether or not the response data type is continuous or categorical
  if (responseDataType == "continuous") {
    
    # Perform naive Bayes on the entire matrix if data is continuous,
    logLikelihoodRatios <- scaledResidualsFilteredToLlr.naiveBayes.evenWidthBins(
      scaledResiduals = scaledResiduals,
      averageSamplesPerBin = averageSamplesPerBin, 
      minFrequencyInTails = minFrequencyInTails)
    
  } else if (responseDataType == "binary" | responseDataType == "ordinal") {
    
    # For every category, perform discretization separately
    # set as binary or ordinal
    
    # Loop through every category, 
    for (categoryValue in unique(actual)) {
      
      # Get a logical vector indicating which rows of the scaled residuals matrix corresponds
      # to the current category value.
      samplesToFitInGaussian <- actual == categoryValue
      
      # Perform a Gaussian naive Bayes method on the rows corresponding 
      # to the current category value.
      logLikelihoodRatios[samplesToFitInGaussian, ] <- scaledResidualsFilteredToLlr.naiveBayes.evenWidthBins(
        scaledResiduals = scaledResiduals, samplesToFitInGaussian = samplesToFitInGaussian,
        averageSamplesPerBin = averageSamplesPerBin, minFrequencyInTails = minFrequencyInTails)
    }
  }
  
  return(logLikelihoodRatios)
}


# Function for converting a matrix of scaled residuals to log likelihood ratios.
# Bins are selected in order to separate the null residuals in equal sized groups.
scaledResidualsToLlr.naiveBayes <- function(scaledResiduals, samplesPerBin = 25) {

  # Extract the null-residuals; 
  # the residuals belonging to the matches that are assumed to be correct.
  nullResiduals <- diag(scaledResiduals)
  
  # Extract the alternative residuals; 
  # the residuals belonging to the matches that are assumed to be sample-swaps.
  alternativeResiduals <- scaledResiduals[lower.tri(scaledResiduals) | upper.tri(scaledResiduals)]
  
  # Store the dimensions of the scaled residual matrix, as well as rownames and colnames.
  returnDimensions <- dim(scaledResiduals)
  returnColumnNames <- colnames(scaledResiduals)
  returnRowNames <- rownames(scaledResiduals)
  
  # Calculate the number of bins given the number of null samples per bin and
  # the total number of null residuals
  nBins <- as.integer(length(nullResiduals) / samplesPerBin)
  
  # Split the null residuals in equal sized tiles, 
  # and get the break points which separate these tiles.
  nullTiles <- ntile(nullResiduals, nBins)
  breaks <- sapply(1:nBins, function(bin) min(nullResiduals[nullTiles == bin]))
  
  # Adapt the breaks to span the entire range of scaled residuals + a buffer of 1.
  breaks[1] <- min(scaledResiduals - 1)
  breaks <- c(breaks, max(scaledResiduals) + 1)
  
  # Split the alternative residuals in tiles according to the same break points
  # that separate the null residuals in equal sized tiles.
  alternativeTiles <- cut(alternativeResiduals, breaks = breaks, labels = FALSE)
  
  # Remove the null and alternative residuals to clear memory.
  rm(alternativeResiduals)
  gc()
  
  # Get the density / likelihood of the null residuals for each of the bins.
  nullLikelihoods <- sapply(
    1:nBins, 
    function(bin) sum(nullTiles == bin) / length(nullTiles))
  
  # Get the density / likelihood of the alternative residuals for each of the bins.
  alternativeLikelihoods <- sapply(
    1:nBins, 
    function(bin) sum(alternativeTiles == bin) / length(alternativeTiles))
  
  # Remove the null and alternative tiles to clear memory.
  rm(nullTiles)
  rm(alternativeTiles)
  gc()
  
  # Determine, for each of the bins, 
  # the ratio between the densities / likelihoods of the alternative compared to the null residuals.
  likelihoodRatioMap <- alternativeLikelihoods / nullLikelihoods
  
  # Apply breaks again on all the residuals.
  allTiles <- cut(scaledResiduals, breaks = breaks, labels = FALSE)
  
  # Remove the scaled residuals from memory.
  rm(scaledResiduals)
  gc()
  
  # Retrieve, for every of the bins / tiles, its respective likelihood ratio.
  likelihoodRatios <- likelihoodRatioMap[allTiles]
  
  # Apply the same dimensions, row names and column names as where used for the input matrix.
  dim(likelihoodRatios) <- returnDimensions
  rownames(likelihoodRatios) <- returnRowNames
  colnames(likelihoodRatios) <- returnColumnNames
  
  # Return the log-transformed likelihood ratios.
  return(log(likelihoodRatios))
}

# Function for converting a selection of a matrix of scaled residuals to likelihood ratios.
# This function employs a Gaussian naive Bayes method to calculate likelihoods.
# Known problems:
# 1. For very low log likelihoods, the value is expected to be -Inf (< -30 approximately)
# 2. If the alternative log likelihood is -Inf as well as the null log likelihood, NaN is returned.
#    When only the null log likelihood is -Inf. Positive infinity is returned.
scaledResidualsFilteredToLlr.gaussianNaiveBayes <- function(scaledResiduals, 
                                                            samplesToFitInGaussian = T,
                                                            error = 0.01) {
  
  # Force the residuals into a normal distribution
  scaledResiduals[lower.tri(scaledResiduals, diag = TRUE) 
                  & upper.tri(scaledResiduals, diag = TRUE)
                  & samplesToFitInGaussian] <- forceNormal(
                    scaledResiduals[lower.tri(scaledResiduals, diag = TRUE) 
                                    & upper.tri(scaledResiduals, diag = TRUE)
                                    & samplesToFitInGaussian])
  
  # Force the residuals into a normal normal distribution
  scaledResiduals[(lower.tri(scaledResiduals) | upper.tri(scaledResiduals)) 
                  & samplesToFitInGaussian] <- forceNormal(
                    scaledResiduals[(lower.tri(scaledResiduals) | upper.tri(scaledResiduals)) 
                                    & samplesToFitInGaussian])
  
  # Extract the null-residuals; 
  # the residuals belonging to the matches that are assumed to be correct.
  nullResiduals <- scaledResiduals[lower.tri(scaledResiduals, diag = TRUE) 
                                   & upper.tri(scaledResiduals, diag = TRUE)
                                   & samplesToFitInGaussian]
  
  # Extract the alternative residuals; 
  # the residuals belonging to the matches that are assumed to be sample-swaps.
  alternativeResiduals <- scaledResiduals[(lower.tri(scaledResiduals) | upper.tri(scaledResiduals)) 
                                          & samplesToFitInGaussian]
  
  # Calculate the mean of both groups of residuals.
  nullMean <- mean(nullResiduals)
  alternativeMean <- mean(alternativeResiduals)
  
  # Calculate the standard deviation in both groups of residuals.
  nullSd <- sd(nullResiduals)
  alternativeSd <- sd(alternativeResiduals)
  
  # Declare function for calculating the difference between two log-transformed
  # likelihoods. returns -Inf if l1 and l2 are (approximately) equal. l1 must be > l2.
  # adapted from: https://stats.stackexchange.com/a/383524 on 22-07-2020.
  # using log(-expm1(...)) in comparison to log1p(-exp(...)) makes the function usable for
  # smaller values.
  subtractLogTransformedProbabilities <- function(l1, l2) {
    # return(l1 + log1p(-exp(-(l1 - l2))))
    return(l1 + log(-expm1(-(l1 - l2))))
  }
  
  # Declare function for calculating the log likelihood of 'value' being sampled from 
  # a normal distribution with the given mean and standard deviation 'sd'.
  gaussianLogLikelihood <- function(value, mean, sd) {
    logLikelihoods <- subtractLogTransformedProbabilities(
      pnorm(value + (error / 2), mean = mean, sd = sd, log.p = TRUE),
      pnorm(value - (error / 2), mean = mean, sd = sd, log.p = TRUE)
    )
    
    # log likelihood is -Inf if log(p) of 
    return(logLikelihoods)
  }
  
  # Calculate, for every residual, the likelihood of the residual being sampled from
  # the corresponding normal distribution of null residuals.
  nullLogLikelihoods <- gaussianLogLikelihood(
    scaledResiduals[samplesToFitInGaussian, ], nullMean, nullSd)
  
  # Calculate, for every residual, the likelihood of the residual being sampled from
  # the corresponding normal distribution of alternative residuals.
  alternativeLogLikelihoods <- gaussianLogLikelihood(
    scaledResiduals[samplesToFitInGaussian, ], alternativeMean, alternativeSd)

  # For every residual, calculate the likelihood ratio the residual belonging to a
  # sample swap.
  logLikelihoodRatios <- alternativeLogLikelihoods - nullLogLikelihoods
  
  return(logLikelihoodRatios)
}

# Function for converting a matrix of scaled residuals to likelihood ratios.
# If the response data type of a trait is categorical or dichotomous, for
# outcome separate Gaussian curves are fitted.
# This function employs a Gaussian naive Bayes method to calculate likelihoods.
# Known problems:
# 1. For very low log likelihoods, the value is expected to be -Inf (< -30 approximately)
# 2. If the alternative log likelihood is -Inf as well as the null log likelihood, NaN is returned.
#    When only the null log likelihood is -Inf. Positive infinity is returned.
scaledResidualsToLlr.gaussianNaiveBayes <- function(
  scaledResiduals, actual, responseDataType = "continuous") {
  
  # Store the dimensions of the scaled residual matrix, as well as rownames and colnames.
  logLikelihoodRatios <- matrix(nrow = nrow(scaledResiduals), 
                                ncol = ncol(scaledResiduals),
                                dimnames = list(rownames(scaledResiduals), colnames(scaledResiduals)))
  
  # Check whether or not the response data type is continuous or categorical
  if (responseDataType == "continuous") {
    
    # Perform Gaussian naive Bayes on the entire matrix if data is continuous,
    # and we thus assume the response data to be normally distributed.
    logLikelihoodRatios <- scaledResidualsFilteredToLlr.gaussianNaiveBayes(scaledResiduals)
    
  } else if (responseDataType == "binary" | responseDataType == "ordinal") {
    
    # For every category, fit separate Gaussian curves if the response data type is either
    # set as binary or ordinal

    # Loop through every category, 
    for (categoryValue in unique(actual)) {
      
      # Get a logical vector indicating which rows of the scaled residuals matrix corresponds
      # to the current category value.
      samplesToFitInGaussian <- actual == categoryValue

      # Perform a Gaussian naive Bayes method on the rows corresponding 
      # to the current category value.
      logLikelihoodRatios[samplesToFitInGaussian, ] <- scaledResidualsFilteredToLlr.gaussianNaiveBayes(
        scaledResiduals = scaledResiduals, samplesToFitInGaussian = samplesToFitInGaussian)
    }
  }
  
  return(logLikelihoodRatios)
}

# Define function for calculating the AUC
calculate.auc <- function(actual, predictor) {
  pr <- prediction(predictor, actual, label.ordering = levels(actual))
  
  auc <- performance(pr, measure = "auc")
  auc <- auc@y.values[[1]]
  return(auc)
}

# Function for forcing normal distribution
forceNormal <- function(x) {

  ranked <- rank(x)
  
  pValues <- (0.5 + ranked - 1.0) / length(ranked)
  
  return(qnorm(pValues, 
               mean = mean(x), 
               sd = sd(x)))
}

# Function for correcting 
# correct <- function(df) {
#   colnames(df) <- c("y", colnames(df)[2:ncol(df)])
#   model.2 <- lm(y ~ . + .^2, df)
#   print(summary(model.2))
#   return(resid(model.2))
# }


plotResiduals <- function(residualsDataFrame, phenotypeTable, responseDataType) {
  
  # Check whether or not the response data type is continuous or categorical
  if (responseDataType == "continuous") {
    
    # Perform Gaussian naive Bayes on the entire matrix if data is continuous,
    # and we thus assume the response data to be normally distributed.
    print(ggplot(residualsDataFrame, aes(x=residuals, stat(density), fill=group)) +
            geom_histogram(bins = 36, alpha=.5, position="identity") +
            geom_rug(data = residualsDataFrame %>% filter(geno != original & group == "null"), 
                     aes(x=residuals), inherit.aes=F) +
            xlab("Unscaled residuals") + ggtitle(paste0("Unscaled residuals for trait '", trait, "'")))
    
  } else if (responseDataType == "binary" | responseDataType == "ordinal") {
    
    # For every category, fit separate Gaussian curves if the response data type is either
    # set as binary or ordinal
    
    # Loop through every category, 
    for (categoryValue in unique(phenotypeTable$VALUE)) {
      
      # Get a logical vector indicating which rows of the scaled residuals matrix corresponds
      # to the current category value.
      samplesToPlot <- phenotypeTable %>%
        filter(VALUE == categoryValue) %>%
        pull(pheno)
      
      residualsToPlot <- residualsDataFrame %>%
        filter(Var1 %in% samplesToPlot)
      
      # Perform a Gaussian naive Bayes method on the rows corresponding 
      # to the current category value.
      print(ggplot(residualsToPlot, aes(x=residuals, stat(density), fill=group)) +
              geom_histogram(bins = 36, alpha=.5, position="identity") +
              geom_rug(data = residualsToPlot %>% filter(geno != original & group == "null"), 
                       aes(x=residuals), inherit.aes=F) +
              xlab("Unscaled residuals") + ggtitle(paste0("Unscaled residuals for trait '", trait, "'")))
    }
  }
}

##############################
# Run
##############################
args <- parser$parse_args(commandArgs(trailingOnly = TRUE))
# args <- parser$parse_args(c("--trait-gwas-mapping", "/groups/umcg-lld/tmp01/other-users/umcg-rwarmerdam/pgs_based_mixup_correction/scripts/r-scripts/pgs_based_sample_mix-up_correction/trait-gwas-mapping.txt",
#                             "--sample-coupling-file", "/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/data/lifelines/processed/pgs.sample-coupling-file.ugli.perm_5120samples_26mixUps.txt",
#                             "--base-pgs-path", "/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/output/PRScs/20200811/",
#                             "--phenotypes-file", "/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/data/lifelines/processed/pgs.phenotypes.ugli.dat",
#                             "--out", "/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/output/sample-swap-prediction/20200811.test/",
#                             "--llr-bayes-method", "efi-discretization", "30"))

message(strwrap(prefix = " ", initial = "", paste(
  "Loading trait-gwas-mapping:\n", args$trait_gwas_mapping)))

# Load table containing paths for the plink output 
# and corresponding phenotype labels.
traitGwasMappingPath <- args$trait_gwas_mapping
traitDescriptionsTable <- fread(
  traitGwasMappingPath, 
  quote="", header=T, sep = "\t",
  col.names=c("trait", "traitDataType", "summaryStatistics", 
              "sampleSizeOfGwas", "numberOfCategories"), 
  stringsAsFactors=F)

message(strwrap(prefix = " ", initial = "", paste(
  "Loading polygenic scores from:\n", args$trait_gwas_mapping)))

# Get the paths to the polygenic scores.
basePathWithPolygenicScores <- args$base_pgs_path
traitDescriptionsTable$polygenicScoreFilePath <- file.path(basePathWithPolygenicScores, 
                                                 traitDescriptionsTable$summaryStatistics, 
                                                 "full.UGLI.pgs.profile")

# Get the output path
out <- args$out

# Set the number of bins to use for the Naive Bayes method.
samplesPerNaiveBayesBin <- 25

# naive bayes method
naiveBayesMethod <- args$llr_bayes_method[1]
stopifnot("method must be one of the following: <ewi-discretization | efi-discretization | gaussian>" = 
            naiveBayesMethod %in% c("ewi-discretization", "efi-discretization", "gaussian"))

if ((naiveBayesMethod == "ewi-discretization" | naiveBayesMethod == "efi-discretization") 
    & length(args$llr_bayes_method) == 2) {
  
  samplesPerNaiveBayesBin <- as.numeric(args$llr_bayes_method[2])
}

message(strwrap(prefix = " ", initial = "", paste(
  "Loading phenotype table:\n", args$phenotypes_file)))

# Load the phenotypes 
phenotypesFilePath <- args$phenotypes_file
phenotypesTable <- fread(phenotypesFilePath, header=T, quote="", sep="\t",
                         col.names = c("ID", "AGE", "SEX", "VALUE", "TRAIT")) %>%
  mutate(SEX = factor(SEX, levels = c("Female", "Male"))) %>%
  group_by(ID) %>%
  filter(!any(AGE < 18)) %>%
  ungroup()

link <- data.frame(geno = unique(phenotypesTable$ID), pheno = unique(phenotypesTable$ID))

# Get the link path
if (!is.null(args$sample_coupling_file)) {
  message(strwrap(prefix = " ", initial = "", paste(
    "Loading sample coupling table:\n", args$sample_coupling_file)))
  
  sampleCouplingFilePath <- args$sample_coupling_file
  link <- fread(sampleCouplingFilePath, stringsAsFactors=F, header=T) %>%
    filter(pheno %in% unique(phenotypesTable$ID))
} else {
  message(strwrap(prefix = " ", initial = "", 
                  "Not using special sample coupling table"))
}

if (FALSE) {
  link <- link %>%
    slice_sample(n = 4096, order_by = geno)
}

aggregatedLlrMatrix <- matrix(nrow = nrow(link), 
                              ncol = nrow(link), 
                              dimnames = list(link$pheno, link$geno), 0)

aggregatedNumberOfTraits <- matrix(nrow = nrow(link), 
                                   ncol = nrow(link), 
                                   dimnames = list(link$pheno, link$geno), 0)

pearson.correlations <- data.frame(trait = traitDescriptionsTable$trait,
                                   pearson.not_corrected = 0.0, 
                                   pearson.corrected.sex = 0.0,
                                   pearson.corrected.age = 0.0,
                                   pearson.corrected.both = 0.0)
rownames(pearson.correlations) <- pearson.correlations$phenotype

# Loop trough traits

for (traitIndex in 1:nrow(traitDescriptionsTable)) {
  
  polygenicScoreFilePath <- traitDescriptionsTable$polygenicScoreFilePath[traitIndex]
  trait <- traitDescriptionsTable$trait[traitIndex]
  responseDataType <- traitDescriptionsTable$traitDataType[traitIndex]
  
  if (!dir.create(file.path(out, trait), recursive = T)) {
    warning(paste0("Could not create directory '", file.path(out, trait), "'"))
  }
  
  phenotypeTable <- phenotypesTable %>%
    filter(TRAIT == trait) %>%
    rename(pheno = ID) %>%
    inner_join(link, by="pheno")
  
  # Give status update
  message(paste0(traitIndex, " / ", nrow(traitDescriptionsTable), 
                 ": '", trait, "' (", responseDataType, ")."))
  
  if (responseDataType == "binary") {
    phenotypeFrequencyTable <- table(phenotypeTable$VALUE)

    message(paste0("    Available for ", nrow(phenotypeTable), 
                   " samples (number of 0's = ", phenotypeFrequencyTable["0"],
                   ", 1's = ", phenotypeFrequencyTable["1"], ")."))
    
    if (any(phenotypeFrequencyTable < 50)) {
      message(paste0("Not enough samples present in group '", names(phenotypeFrequencyTable)[phenotypeFrequencyTable < 50], "'. Skipping..."))
      next
    }
    
  } else if (responseDataType == "continuous") {
    message(paste0("    Available for ", nrow(phenotypeTable), " samples."))
  }

  message(paste0("    Loading polygenic scores from '", polygenicScoreFilePath, "'..."))
  
  # Read the PLINK polygenic score table.
  polygenicScores <- read.table(
    polygenicScoreFilePath,
    header=T) %>%
    rename(PGS = SCORESUM)

  completeTable <- phenotypeTable %>%
    inner_join(polygenicScores, by = c("geno" = "IID"))
  
  initial.cor.test.results <- cor.test(completeTable$VALUE, completeTable$PGS)
  message(paste0("Initial R-squared of correlation = ", initial.cor.test.results$estimate ^ 2))
  pearson.correlations[traitIndex, "pearson.not_corrected"] <- initial.cor.test.results$estimate
  
  # Correct for age and sex
  trait2pgs.corrected <- completeTable
  model.complete <- lm(VALUE ~ AGE + SEX + AGE * SEX, data = trait2pgs.corrected)
  trait2pgs.corrected$actual <- resid(model.complete)
  
  # Output the correlation of the corrected traits
  corrected.cor.test.results <- cor.test(trait2pgs.corrected$actual, trait2pgs.corrected$PGS)
  message(paste0("R-squared of corrected traits = ", corrected.cor.test.results$estimate ^ 2))
  pearson.correlations[traitIndex, "pearson.corrected.both"] <- corrected.cor.test.results$estimate

  # Calculate z-score matrix based on polygenic scores and actual phenotypes,
  # using the chosen function for calculating residuals.
  
  pdf(file.path(out, trait, "/debugFigures.pdf"))
  par(xpd = NA)
  
  # Calculate the scaled residuals for every combination
  scaledResidualsMatrix <- calculate.scaledResiduals(
    estimate = completeTable$PGS, 
    actual = completeTable$VALUE, 
    covariates = completeTable[c("AGE", "SEX")],
    responseDataType = responseDataType)
  
  dev.off()
  
  # The rows correspond to phenotype samples, the cols to genotype samples (polygenic scores)
  rownames(scaledResidualsMatrix) <- completeTable$pheno
  colnames(scaledResidualsMatrix) <- completeTable$geno
  
  # Write scaled residuals matrix
  write.table(scaledResidualsMatrix, file.path(out, trait, "/scaledResidualMatrix.tsv"), 
              sep = "\t", col.names = T, row.names = T, quote = F)

  message("    completed calculating scaled residuals")
  
  message("    calculating log likelihood ratios")
  
  logLikelihoodRatios <- NULL
  
  # Calculate the likelihood ratios for every residual being from the distribution of possible mix-ups.
  if (naiveBayesMethod == "gaussian") {
    
    logLikelihoodRatios <- scaledResidualsToLlr.gaussianNaiveBayes(
      scaledResiduals = scaledResidualsMatrix, actual = completeTable$VALUE, responseDataType = responseDataType)
    
  } else if (naiveBayesMethod == "ewi-discretization") {
    
    logLikelihoodRatios <- scaledResidualsToLlr.naiveBayes.evenWidthBins(
      scaledResiduals = scaledResidualsMatrix, actual = completeTable$VALUE, 
      responseDataType = responseDataType, averageSamplesPerBin = samplesPerNaiveBayesBin)
    
  } else if (naiveBayesMethod == "efi-discretization") {
    
    logLikelihoodRatios <- scaledResidualsToLlr.naiveBayes(
      scaledResiduals = scaledResidualsMatrix,
      samplesPerBin = samplesPerNaiveBayesBin)
    
  }
  
  # Resolve log likelihood ratios that are NaN (not a number)
  llrIsNan <- is.nan(logLikelihoodRatios)
  logLikelihoodRatios[llrIsNan] <- 0
  
  aggregatedNumberOfTraits[rownames(logLikelihoodRatios), colnames(logLikelihoodRatios)] <- 
    aggregatedNumberOfTraits[rownames(logLikelihoodRatios), colnames(logLikelihoodRatios)] + !llrIsNan
  
  write.table(logLikelihoodRatios, file.path(out, trait, "/logLikelihoodRatios.tsv"), 
              sep = "\t", col.names = T, row.names = T, quote = F)

  message("    completed calculating log likelihood ratios.")
  message(paste0("Calculated AUC: ", calculate.auc(
    logLikelihoodRatios$group, lrProducts$scaledLlr)))
  message("    summing log likelihood ratios with aggregated log likelihood ratios...")
  
  aggregatedLlrMatrix[rownames(logLikelihoodRatios), colnames(logLikelihoodRatios)] <- 
    aggregatedLlrMatrix[rownames(logLikelihoodRatios), colnames(logLikelihoodRatios)] + logLikelihoodRatios
  
  message("    aggregation done!")
  
  scaledResidualsDataFrame <- 
    as.data.frame.table(scaledResidualsMatrix, responseName = "scaledResiduals") %>%
    inner_join(link, by = c("Var1" = "pheno")) %>%
    mutate(group = case_when(Var2 == geno ~ "null",
                             TRUE ~ "alternative"))
  
  scaledResidualsDataFrame$group <- factor(scaledResidualsDataFrame$group, levels = c("null", "alternative"))
  
  rm(scaledResidualsMatrix)
  gc()
  
  pdf(file.path(out, trait, "unscaledResidualsHistogram.pdf"),
      width=8, height = 6, useDingbats = FALSE)

  par(xpd = NA)

  plotResiduals(scaledResidualsDataFrame, phenotypeTable, responseDataType)

  dev.off()
  
  # ggplot(scaledResidualsDataFrame, aes(x=scaledResiduals, stat(density), fill=group)) +
  #   geom_histogram(bins = 72, alpha=.5, position="identity") +
  #   geom_rug(data = scaledResidualsDataFrame %>% filter(geno != original & group == "null"), aes(x=scaledResiduals), inherit.aes=F) +
  #   xlab("Unscaled residuals") + ggtitle(paste0("Unscaled residuals for trait '", trait, "'"))
  
  # ggsave("unscaledResidualsHistogram.png", width=8, height=7)
  
  rm(logLikelihoodRatios)
  gc()
}

write.table(aggregatedLlrMatrix, file.path(out, "/aggregatedLogLikelihoodRatiosMatrix.tsv"),
            sep="\t", col.names = T, row.names = T, quote = F)

write.table(aggregatedNumberOfTraits, file.path(out, "/aggregatedNumberOfTraitsMatrix.tsv"),
            sep="\t", col.names = T, row.names = T, quote = F)

lrProducts <- 
  as.data.frame.table(aggregatedLlrMatrix, responseName = "logLikelihoodRatios") %>%
  inner_join(link, by = c("Var1" = "pheno"))

write.table(lrProducts, file.path(out, "/aggregatedLogLikelihoodRatiosDataFrame.tsv"),
            sep="\t", col.names = T, row.names = F, quote = F)

lrProducts$group <- "alternative"
lrProducts$group[lrProducts$original == lrProducts$Var2] <- "null"
lrProducts$group <- factor(lrProducts$group, levels = c("null", "alternative"))

lrProducts <- lrProducts %>% 
  group_by(Var1) %>%
  mutate(scaledLlr = scale(logLikelihoodRatios))

message(paste0("Calculated overall AUC: ", calculate.auc(
  lrProducts$group, lrProducts$scaledLlr)))

permutationTestDataFrame <- lrProducts %>%
  filter(geno == Var2)

message(paste0("Confined AUC: ", calculate.auc(
  permutationTestDataFrame$group, 
  permutationTestDataFrame$scaledLlr)))

ggplot(lrProducts, aes(x=logLikelihoodRatios, stat(density), fill=group)) +
  geom_histogram(bins = 32, alpha=.5, position="identity") +
  xlab("Log likelihood ratios") + ggtitle(paste0("LR overall"))

ggsave(file.path(out, "/likelihoodRatioHistogram.png"), width=8, height=7)

ggplot(lrProducts, aes(x=group, y=logLikelihoodRatios)) +
  geom_boxplot() + ggtitle(paste0("Log likelihood ratio distributions comparison overall"))

ggsave(file.path(out, "/likelihoodRatioBoxplot.png"), width=8, height=7)

numberOfTraits <- 
  as.data.frame.table(aggregatedNumberOfTraits, responseName = "numberOfTraits")

phenoSamples <- unique(lrProducts[lrProducts$geno != lrProducts$original, "Var1"])

# phenoSamples <- sample(unique(lrProducts[lrProducts$geno == lrProducts$original, "Var1"]), size = 26)
# 
# sampledLrProducts <- lrProducts %>% 
#   #inner_join(numberOfTraits, by=c("Var1", "Var2")) %>%
#   mutate(colourGroup = case_when(
#     Var2 == geno ~ 1,
#     Var2 == original ~ 2,
#     TRUE ~ 3)) %>%
#   filter(Var1 %in% phenoSamples)
# 
# cols = c("2" = "blue", "1" = "red", "3" = "black")
# 
# ggplot(sampledLrProducts %>% filter(is.finite(logLikelihoodRatios)), 
#        aes(y = Var1)) +
#   
#   # Add ridge lines
#   geom_density_ridges(
#     aes(x = logLikelihoodRatios, 
#         point_alpha = as.numeric(colourGroup != 3),
#         point_color = colourGroup),
#     jittered_points = TRUE,
#     position = position_points_jitter(width = 0, height = 0),
#     point_shape = '|', point_size = 3, alpha = 0.7) +
#   scale_point_color_continuous(low = "#0072B2", high = "#D55E00") +
#   #scale_discrete_manual(values = cols, aesthetics = "point_color") +
#   #scale_discrete_manual("point_color", values = cols) +
# 
#   # Add annotation with number of traits, and the number of filtered likelihood ratios
#   # geom_text(
#   #   data=sampledLrProducts %>% group_by(Var1) %>%
#   #     summarise(numberOfTraits = median(numberOfTraits),
#   #               numberOfFiltered = sum(!is.finite(logLikelihoodRatios))),
#   #   position=position_nudge(y=0.64), colour="red", size=3.5,
#   #   hjust = "inward", x = 0,
#   #   aes(label = sprintf("n traits: %d, n filtered: %d", numberOfTraits, numberOfFiltered))) +
# 
#   # Set theme
#   theme_minimal() +
#   theme(legend.position = "none")
# 
# ggsave("~/pgs_based_mixup_correction-ugli/output/sample-swap-prediction/20200811.20201012.efi-discretization.30/ridges.png", width=8, height=20)

# Pearson correlations
# pearson.correlations <- pearson.correlations[order(pearson.correlations$pearson.corrected),]
# pearson.correlations$r2.not_corrected <- pearson.correlations$pearson.not_corrected ^ 2
# pearson.correlations$r2.corrected.sex <- pearson.correlations$pearson.corrected.sex ^ 2
# pearson.correlations$r2.corrected.age <- pearson.correlations$pearson.corrected.age ^ 2
# pearson.correlations$r2.corrected.both <- pearson.correlations$pearson.corrected.both ^ 2
# 
# write.table(pearson.correlations, file.path(out, "/correlations_comparison.tsv"))
# 
# pearson.correlations.melted <- melt(pearson.correlations[,c('phenotype','r2.not_corrected','r2.corrected.sex', 'r2.corrected.age', 'r2.corrected.both')], value.name = "r2")
# pearson.correlations.melted <- pearson.correlations.melted[order(pearson.correlations.melted[,"r2"]),]
# 
# ggplot(pearson.correlations.melted,aes(x = phenotype,y = r2)) +
#   geom_bar(aes(fill = variable),stat = "identity",position = "dodge") +
#   coord_flip()
# 
# ggsave(file.path(out, "/correlations_comparison.png"), width=8, height=7)

# # Calculate the correlation of residual Z-scores between phenotypes
# res.zscore.correlations <- sapply(lRdataFrameList, function(x) sapply(lRdataFrameList, function(y) cor(x["likelihoodRatios"], y["likelihoodRatios"])))
# rownames(res.zscore.correlations) <- phenotypes
# colnames(res.zscore.correlations) <- phenotypes
# res.zscore.correlations <- as.matrix(res.zscore.correlations)
# 
# res.zscore.correlations.melted <- melt(res.zscore.correlations, value.name="pearson.cor")
# 
# # Color Brewer palette
# ggplot(data = res.zscore.correlations.melted, aes(x=Var1, Var2, fill=pearson.cor)) +
#   ggtitle("Correlations of LRs between phenotypes") +
#   xlab("Phenotypes") +
#   ylab("Phenotypes") +
#   geom_tile() +
#   scale_fill_distiller(palette = "Blues") +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))
# 
# ggsave(paste0(dirname(filename), "/residual_LikelihoodRatioCorrelations.png"), width=8, height=7)

# sampleTraitsSummary <- sampledLrProducts %>% group_by(Var1) %>%
#       summarise(numberOfTraits = median(numberOfTraits),
#                 numberOfFiltered = sum(!is.finite(logLikelihoodRatios)))
# 
# lrProducts_Prev <- 
#   as.data.frame.table(aggregatedLlrMatrix, responseName = "logLikelihoodRatios") %>%
#   inner_join(link, by = c("Var1" = "pheno"))
# 
# lrProducts_Prev$group <- "alternative"
# lrProducts_Prev$group[lrProducts_Prev$original == lrProducts_Prev$Var2] <- "null"
# lrProducts_Prev$group <- factor(lrProducts_Prev$group, c("alternative", "null"))
