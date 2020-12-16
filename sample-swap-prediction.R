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
library(pROC)

##############################
# Define argument parser
##############################

parser <- ArgumentParser(description='')
parser$add_argument('--debug', action='store_true', 
                    dest="debug", help="write intermediate llr and residuals files.")
parser$add_argument('--base-fit-model-path',
                    help='path to a directory to write fitted model parmaters to, or to load fitted models from.')
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
  # predictedProbabilities[
  #   cbind(seq_len(length(actual)), actual)] <- 0
  
  # Get the column for which the probability is largest,
  # excluding the column corresponding to the actual class
  # competingProbability <- max.col(predictedProbabilities, "first")
  # Get the corresponding classes as an ordered factor
  # competingClass <- factor(levels(actual)[competingProbability], 
  #                          levels = levels(actual), ordered = T)

  # Calculate the deviance residuals from the probability for the actual class
  deviance <- sqrt(2 * -log(predictedProbabilityForActualClass))
  # Add direction:
  # When the probabilities gravitate towards > actual class: negative direction
  # When the probabilities gravitate towards < actual class: positive direction
  # deviance <- ((actual > competingClass) - (competingClass > actual)) * deviance
  
  return(deviance)
}

# Function that returns a function for calculating residuals.
# The model determining residuals is dependent on the type of response data type.
residualsFunConstructor <- function(estimate, actual, covariates, 
                                    responseDataType = "continuous", useOrderedLogit = TRUE,
                                    modelPath = NULL) {
  
  # If the path to write/read fitted models to already points to an existing file, load this.
  # In this case residual calculation should be performed using the loaded model
  if (!is.null(modelPath) && 0 != length(modelPath) && file.exists(modelPath)) {
    load(modelPath)
  }
  
  # Return a linear model in case 'actual', is a continuous data type.
  if (responseDataType == "continuous" | (responseDataType == "ordinal" & !useOrderedLogit) 
      | exists("olsModel")) {
    
    if (!exists("olsModel")) {
      olsModel <- lm(actual ~ estimate + . + .^2, data = covariates)
      if (0 != length(modelPath)) {
        save(olsModel, file = modelPath)
      }
    }
    
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
  } else if (responseDataType == "binary"
             | exists("logitModel")) {
    
    if (!exists("logitModel")) {
      logitModel <- glm(actual ~ estimate + . + .^2, 
                        family = binomial(link='logit'),
                        data = covariates)
      if (0 != length(modelPath)) {
        save(logitModel, file = modelPath)
      }
    }
    
    print(summary(logitModel))
    
    #plotSigmoid(estimate = estimate, actual = actual, covariates = covariates, logitModel = logitModel)
    
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
  } else if ((responseDataType == "ordinal" & useOrderedLogit)
             | exists("orderedLogitModel")) {
    
    classes <- sort(unique(actual))
    actualOrdered <- factor(actual, levels = classes, ordered = TRUE)
    
    if (!exists("orderedLogitModel")) {
      orderedLogitModel <- polr(actualOrdered ~ estimate + . + .^2, 
                                method = "logistic", Hess = TRUE,
                                data = covariates)
      if (0 != length(modelPath)) {
        save(orderedLogitModel, file = modelPath)
      }
    }
    
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
calculate.scaledResiduals <- function(estimate, actual, covariates, responseDataType, modelPath, useOrderedLogit = TRUE) {
  # Get a function that returns the deviance of an observation to the regression line.
  # This will either be the 'deviance residuals' from a logistic model when the response data type is binary,
  # or this will be the regular residuals corresponding to a linear model.
  residualsFun <- residualsFunConstructor(estimate = estimate,
                                          actual = actual,
                                          covariates = covariates,
                                          responseDataType = responseDataType, 
                                          useOrderedLogit = useOrderedLogit,
                                          modelPath = modelPath)

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
  
  if (nBins == 2) {
    
    breaks <- c(min(x) - 1, mean(c(min(x), max(x))), max(x) + 1)
    
  } else if (n >= (minFrequencyInTails * 2)) {
    
    upperTailLowerBound <- sort(x, partial = n - minFrequencyInTails)[n - minFrequencyInTails]
    upperTailUpperBound <- max(x) + 1
    
    lowerTailUpperBound <- sort(x, partial = minFrequencyInTails)[minFrequencyInTails]
    lowerTailLowerBound <- min(x) - 1
    
    # Separate the space within the outermost bins into a nBins - 2
    binSize <- (upperTailLowerBound - lowerTailUpperBound) / (nBins - 2)
    breaks <- unique(c(lowerTailLowerBound, lowerTailUpperBound,
                lowerTailUpperBound + binSize * (1:(nBins - 3)), 
                upperTailLowerBound, upperTailUpperBound))
  
    tiles <- cut(x, breaks = breaks)
    
    while (min(table(tiles)) < minFrequencyInTails) {
      
      minTile <- which.min(table(tiles))
      breaks[minTile] <- mean(breaks[c(minTile, minTile + 1)])
      breaks <- unique(breaks[-(minTile + 1)])
      
      tiles <- cut(x, breaks = breaks)
    }
    
  } else {
    
    breaks <- c(min(x) - 1, max(x) + 1)
  }
  
  return(unique(breaks))
}

# Function for converting a matrix of scaled residuals to log likelihood ratios.
# Bins are selected while maintaining a regular bandwidth
fitEwiDiscretizationParameters <- function(
  nullValues, alternativeValues, 
  averageSamplesPerBin = 25, minFrequencyInTails = 10) {
  
  # Calculate the number of bins given the number of null samples per bin and
  # the total number of null residuals
  nBins <- as.integer(round(length(nullValues) / averageSamplesPerBin))

  message("Obtaining breaks")
  
  breaks <- adaptedEqualWidthIntervals(nullValues, nBins, minFrequencyInTails)
  breaks[1] <- -Inf
  breaks[length(breaks)] <- Inf
  
  message(paste0("Breaks created: ", length(breaks)))

  nullTiles <- cut(nullValues, breaks = breaks, labels = FALSE)
  alternativeTiles <- cut(alternativeValues, breaks = breaks, labels = FALSE)
  
  message("Breaks applied")

  # Get the density / likelihood of the null residuals for each of the bins.
  nullLikelihoods <- sapply(
    1:max(nullTiles), 
    function(bin) sum(nullTiles == bin) / length(nullTiles))
  
  # Get the density / likelihood of the alternative residuals for each of the bins.
  alternativeLikelihoods <- sapply(
    1:max(nullTiles), 
    function(bin) sum(alternativeTiles == bin) / length(alternativeTiles))
  
  # Remove the null and alternative tiles to clear memory.
  rm(nullTiles)
  rm(alternativeTiles)
  gc()
  
  # Determine, for each of the bins, 
  # the ratio between the densities / likelihoods of the alternative compared to the null residuals.
  likelihoodRatioMap <- alternativeLikelihoods / nullLikelihoods
  
  print(any(is.infinite(likelihoodRatioMap)))
  
  return(list(breaks = breaks, likelihoodRatioMap = likelihoodRatioMap))
}

# Function for converting a matrix of scaled residuals to log likelihood ratios.
# Bins are selected in order to separate the null residuals in equal sized groups.
fitEfiDiscretizationParameters <- function(
  nullValues, alternativeValues, samplesPerBin = 25) {
  
  # Calculate the number of bins given the number of null samples per bin and
  # the total number of null residuals
  nBins <- as.integer(length(nullValues) / samplesPerBin)
  
  # Split the null residuals in equal sized tiles, 
  # and get the break points which separate these tiles.
  nullTiles <- ntile(nullValues, nBins)
  breaks <- sapply(1:nBins, function(bin) min(nullValues[nullTiles == bin]))
  
  # Adapt the breaks to span the entire range of values + a buffer of 1.
  breaks[1] <- -Inf
  breaks <- c(breaks, Inf)
  
  # Split the alternative residuals in tiles according to the same break points
  # that separate the null residuals in equal sized tiles.
  alternativeTiles <- cut(alternativeValues, breaks = breaks, labels = FALSE)
  
  # Remove the null and alternative values to clear memory.
  rm(nullValues)
  rm(alternativeValues)
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
  
  return(list(breaks = breaks, likelihoodRatioMap = likelihoodRatioMap))
}

calculateDiscritizedLoglikelihoodRatios <- function(
  values, parameters) {
  
  # Apply breaks again on all the residuals.
  allTiles <- cut(values, breaks = parameters$breaks, labels = FALSE)
  
  # Remove the scaled residuals from memory.
  rm(values)
  gc()
  
  # Retrieve, for every of the bins / tiles, its respective likelihood ratio.
  likelihoodRatios <- parameters$likelihoodRatioMap[allTiles]
  
  # Return the log-transformed likelihood ratios.
  return(log(likelihoodRatios))
}

# Function for converting a selection of a matrix of scaled residuals to likelihood ratios.
# This function employs a Gaussian naive Bayes method to calculate likelihoods.
# Known problems:
# 1. For very low log likelihoods, the value is expected to be -Inf (< -30 approximately)
# 2. If the alternative log likelihood is -Inf as well as the null log likelihood, NaN is returned.
#    When only the null log likelihood is -Inf. Positive infinity is returned.
gaussianNaiveBayes <- function(values, 
                               parameters,
                               error = 0.01) {
  
  # Calculate the mean of both groups of residuals.
  nullMean <- parameters$null$mean
  alternativeMean <- parameters$alternative$mean
  
  # Calculate the standard deviation in both groups of residuals.
  nullSd <- parameters$null$sd
  alternativeSd <- parameters$alternative$sd
  
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
  # nullLogLikelihoods <- gaussianLogLikelihood(
  #   values, nullMean, nullSd)
  nullLogLikelihoods <- dnorm(
    values, nullMean, nullSd, log = TRUE)
  
  # Calculate, for every residual, the likelihood of the residual being sampled from
  # the corresponding normal distribution of alternative residuals.
  # alternativeLogLikelihoods <- gaussianLogLikelihood(
  #   values, alternativeMean, alternativeSd)
  alternativeLogLikelihoods <- dnorm(
    values, alternativeMean, alternativeSd, log = TRUE)

  # For every residual, calculate the likelihood ratio the residual belonging to a
  # sample swap.
  logLikelihoodRatios <- alternativeLogLikelihoods - nullLogLikelihoods
  
  return(logLikelihoodRatios)
}

# Function that fits a Gaussian model, and returns the fitted parameters
fitGaussian <- function(values) {
  return(list(mean = mean(values),
       sd = sd(values)))
}

# Function that calculates log likelihood ratios for a selection of a matrix
calculate.logLikelihoodRatiosForSelection <- function(
  valueMatrix, selection = T,
  naiveBayesMethod, samplesPerBin, classifierPath = NULL) {
  
  # Filter value matrix with selection
  valueMatrixFiltered <- valueMatrix[selection, ]
  
  # Store the dimensions of the value matrix, as well as row names and column names.
  returnDimensions <- dim(valueMatrixFiltered)
  returnColumnNames <- colnames(valueMatrixFiltered)
  returnRowNames <- rownames(valueMatrixFiltered)
  
  # Extract the null-residuals; 
  # the residuals belonging to the matches that are assumed to be correct.
  nullValues <- valueMatrix[lower.tri(valueMatrix, diag = TRUE) 
                                   & upper.tri(valueMatrix, diag = TRUE)
                                   & selection]
  
  # Extract the alternative residuals; 
  # the residuals belonging to the matches that are assumed to be sample-swaps.
  alternativeValues <- valueMatrix[(lower.tri(valueMatrix) | upper.tri(valueMatrix)) 
                                          & selection]
  
  rm(valueMatrix)
  gc()
  
  # Store the dimensions of the scaled residual matrix, as well as rownames and colnames.
  logLikelihoodRatios <- matrix(nrow = nrow(valueMatrixFiltered), 
                                ncol = ncol(valueMatrixFiltered),
                                dimnames = list(rownames(valueMatrixFiltered), 
                                                colnames(valueMatrixFiltered)))
  
  # If the path to write/read fitted models to already points to an existing file, load this.
  # In this case residual calculation should be performed using the loaded model
  if (!is.null(classifierPath) 
      && 0 != length(classifierPath) 
      && file.exists(classifierPath)) {
    load(classifierPath)
  }
  
  # Loop through categories
  # Calculate the likelihood ratios for every residual being from the distribution of possible mix-ups.
  if (naiveBayesMethod == "gaussian" | exists("gaussianParameters")) {
    message("Using gaussian naive bayes method")
    
    if (!exists("gaussianParameters")) {
      gaussianParameters <- list(null = fitGaussian(nullValues), 
                            alternative = fitGaussian(alternativeValues))
      if (0 != length(classifierPath)) {
        save(gaussianParameters, file = classifierPath)
        message(paste0("Saved fitted classifier to following path: ", classifierPath))
      }
    } else {
      message("Using prefitted classifier")
    }

    logLikelihoodRatios <- gaussianNaiveBayes(
      valueMatrixFiltered, gaussianParameters)
    rm(gaussianParameters)
    
  } else if (naiveBayesMethod == "ewi-discretization" | exists("ewiDiscretizationParameters")) {
    message("Using naive bayes method with ewi discretization")
    
    if (!exists("ewiDiscretizationParameters")) {
      ewiDiscretizationParameters <- fitEwiDiscretizationParameters(
        nullValues = nullValues, alternativeValues = alternativeValues, 
        averageSamplesPerBin = samplesPerBin)
      if (0 != length(classifierPath)) {
        save(ewiDiscretizationParameters, file = classifierPath)
        message(paste0("Saved fitted classifier to following path: ", classifierPath))
      }
    } else {
      message("Using prefitted classifier")
    }
    
    print(ewiDiscretizationParameters)
    
    logLikelihoodRatios <- calculateDiscritizedLoglikelihoodRatios(
      valueMatrixFiltered, ewiDiscretizationParameters)
    rm(ewiDiscretizationParameters)
    
  } else if (naiveBayesMethod == "efi-discretization" | exists("efiDiscretizationParameters")) {
    message("Using naive bayes method with efi discretization")
    
    if (!exists("efiDiscretizationParameters")) {
      efiDiscretizationParameters <- fitEfiDiscretizationParameters(
        nullValues = nullValues, alternativeValues = alternativeValues, samplesPerBin = samplesPerBin)
      if (0 != length(classifierPath)) {
        save(efiDiscretizationParameters, file = classifierPath)
        message(paste0("Saved fitted classifier to following path: ", classifierPath))
      }
    } else {
      message("Using prefitted classifier")
    }
    
    logLikelihoodRatios <- calculateDiscritizedLoglikelihoodRatios(
      valueMatrixFiltered, efiDiscretizationParameters)
    rm(efiDiscretizationParameters)
    
  }
  
  # Apply the same dimensions, row names and column names as where used for the input matrix.
  dim(logLikelihoodRatios) <- returnDimensions
  rownames(logLikelihoodRatios) <- returnRowNames
  colnames(logLikelihoodRatios) <- returnColumnNames
  
  return(logLikelihoodRatios)
}

# Function that calculates log likelihoods given a naive Bayes method.
calculate.logLikelihoodRatios <- function(
  valueMatrix, actual, responseDataType, naiveBayesMethod, samplesPerBin, classifierPath) {
  
  # Store the dimensions of the scaled residual matrix, as well as rownames and colnames.
  logLikelihoodRatios <- matrix(nrow = nrow(valueMatrix), 
                                ncol = ncol(valueMatrix),
                                dimnames = list(rownames(valueMatrix), colnames(valueMatrix)))
  
  message("Starting calculating log likelihood ratios...")
  
  # Check whether or not the response data type is continuous or categorical
  if (responseDataType == "continuous") {
    
    # Perform naive Bayes on the entire matrix if data is continuous,
    logLikelihoodRatios <- calculate.logLikelihoodRatiosForSelection(
      valueMatrix = valueMatrix,
      naiveBayesMethod = naiveBayesMethod,
      samplesPerBin = samplesPerBin,
      classifierPath = file.path(classifierPath,
        paste0("likelihoodClassifier.rda")))
    
  } else if (responseDataType == "binary" | responseDataType == "ordinal") {
    
    # For every category, perform discretization separately
    # set as binary or ordinal
    
    message(paste0(
      "Performing log likelihood calculation separately for categories: (", 
      paste0(unique(actual), collapse = ", "), ")"))
    
    # Loop through every category, 
    for (categoryValue in unique(actual)) {
      
      # Get a logical vector indicating which rows of the scaled residuals matrix corresponds
      # to the current category value.
      selection <- actual == categoryValue
      
      message(paste0("  (", categoryValue, ", ", sum(selection), " / ", length(selection), ")"))
      
      # Perform a Gaussian naive Bayes method on the rows corresponding 
      # to the current category value.
      logLikelihoodRatios[selection, ] <- calculate.logLikelihoodRatiosForSelection(
        valueMatrix = valueMatrix,
        selection = selection,
        naiveBayesMethod = naiveBayesMethod,
        samplesPerBin = samplesPerBin,
        classifierPath = file.path(classifierPath, 
          paste0("likelihoodClassifier", "_category", categoryValue, ".rda")))
          
      message("")
    }
  }
  
  message("Done!")
  
  return(logLikelihoodRatios)
}

# Function for forcing normal distribution
forceNormal <- function(x) {

  ranked <- rank(x)
  
  pValues <- (0.5 + ranked - 1.0) / length(ranked)
  
  return(qnorm(pValues, 
               mean = mean(x), 
               sd = sd(x)))
}

plotResiduals <- function(residualsDataFrame, phenotypeTable, responseDataType) {
  
  # Check whether or not the response data type is continuous or categorical
  if (responseDataType == "continuous") {
    
    # Perform Gaussian naive Bayes on the entire matrix if data is continuous,
    # and we thus assume the response data to be normally distributed.
    print(ggplot(residualsDataFrame, aes(x=scaledResiduals, stat(density), fill=group)) +
            geom_histogram(bins = 36, alpha=.5, position="identity") +
            geom_rug(data = residualsDataFrame %>% filter(geno != original & group == "null"), 
                     aes(x=scaledResiduals), inherit.aes=F) +
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
      print(ggplot(residualsToPlot, aes(x=scaledResiduals, stat(density), fill=group)) +
              geom_histogram(bins = 36, alpha=.5, position="identity") +
              geom_rug(data = residualsToPlot %>% filter(geno != original & group == "null"), 
                       aes(x=scaledResiduals), inherit.aes=F) +
              xlab("Unscaled residuals") + ggtitle(paste0("Unscaled residuals for trait '", trait, "'", "(cat: ", categoryValue, ")")))
    }
  }
}

##############################
# Run
##############################
# args <- parser$parse_args(c("--trait-gwas-mapping", "/groups/umcg-lld/tmp01/other-users/umcg-rwarmerdam/pgs_based_mixup_correction/scripts/r-scripts/pgs_based_sample_mix-up_correction/trait-gwas-mapping.txt",
#                             "--sample-coupling-file", "/home/umcg-rwarmerdam/pgs_based_mixup_correction-ugli/data/lifelines/processed/pgs.sample-coupling-file.ugli.20201120.10080samples.txt",
#                             "--base-pgs-path", "/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/output/PRScs/20201120/",
#                             "--phenotypes-file", "/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/data/lifelines/processed/pgs.phenotypes.ugli.dat",
#                             "--out", "/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/output/sample-swap-prediction/20200811.test/",
#                             "--llr-bayes-method", "NA", "50",
#                             "--base-fit-model-path", "/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/output/sample-swap-prediction/fitted-paramters/"))
args <- parser$parse_args(commandArgs(trailingOnly = TRUE))

message(strwrap(prefix = " ", initial = "", paste(
  "Loading trait-gwas-mapping:\n", args$trait_gwas_mapping)))

# Load table containing paths for the plink output 
# and corresponding phenotype labels.
traitGwasMappingPath <- args$trait_gwas_mapping
traitDescriptionsTable <- fread(
  traitGwasMappingPath, 
  quote="", header=T, sep = "\t", 
  stringsAsFactors=F) %>%
  filter(trait == "EduYears")

message(strwrap(prefix = " ", initial = "", paste(
  "Loading polygenic scores from:\n", args$trait_gwas_mapping)))

# Get the paths to the polygenic scores.
basePathWithPolygenicScores <- args$base_pgs_path
traitDescriptionsTable$polygenicScoreFilePath <- file.path(basePathWithPolygenicScores, 
                                                 traitDescriptionsTable$summaryStatistics, 
                                                 "full.UGLI.pgs.profile")

# Get the output path
out <- args$out

debug <- args$debug

# Set the number of bins to use for the Naive Bayes method.
samplesPerNaiveBayesBin <- 25

# Set the likelihood alpha
likelihoodRatioDifferenceAlpha <- 0.05

# naive bayes method
naiveBayesMethod <- args$llr_bayes_method[1]
stopifnot("method must be one of the following: <ewi-discretization | efi-discretization | gaussian>" = 
            naiveBayesMethod %in% c("ewi-discretization", "efi-discretization", "gaussian", "NA"))

if ((naiveBayesMethod == "ewi-discretization" | naiveBayesMethod == "efi-discretization" | naiveBayesMethod == "NA") 
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

link <- data.frame(geno = unique(phenotypesTable$ID), pheno = unique(phenotypesTable$ID), stringsAsFactors = F)

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

predictingInducedMixUps <- F

if (!("original" %in% colnames(link))) {
  link$original <- link$geno
} else {
  predictingInducedMixUps <- T
}

modelBasePath <- NULL

if (!is.null(args$base_fit_model_path)) {
  modelBasePath <- args$base_fit_model_path
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

traitDescriptionsTable <- traitDescriptionsTable %>%
  mutate(confinedAuc = NA_real_,
         matrixWideAuc = NA_real_,
         pValue = NA_real_,
         traitOutputDir = NA_character_,
         modelBasePath = modelBasePath) %>%
  group_by(trait) %>%
  mutate(naiveBayesMethod = case_when(
    !is.na(naiveBayesMethod) 
    & naiveBayesMethod %in% c("gaussian, efi-discretization", "ewi-discretization") ~ naiveBayesMethod,
    traitDataType == "continuous" ~ "gaussian",
    traitDataType %in% c("ordinal", "binary") ~ "efi-discretization"),
    samplesPerNaiveBayesBin = samplesPerNaiveBayesBin)
  
# Loop trough traits

for (traitIndex in 1:nrow(traitDescriptionsTable)) {
  
  polygenicScoreFilePath <- traitDescriptionsTable$polygenicScoreFilePath[traitIndex]
  trait <- traitDescriptionsTable$trait[traitIndex]
  responseDataType <- traitDescriptionsTable$traitDataType[traitIndex]
  naiveBayesMethod <- traitDescriptionsTable$naiveBayesMethod[traitIndex]
  
  traitFileName <- paste(traitIndex, gsub(" ", "_", trait), sep = ".")
  traitDescriptionsTable[traitIndex, "traitOutputDir"] <- traitFileName
  
  modelPath <- NULL
  pgsPhenotypeModel <- NULL
  
  if (!is.null(modelBasePath)) {
    modelPath <- file.path(modelBasePath, traitFileName)
    pgsPhenotypeModel <- file.path(modelPath, "pgsPhenotypeModel.rda")
  }

  if (!dir.create(file.path(out, traitFileName), recursive = T)) {
    warning(paste0("Could not create directory '", file.path(out, traitFileName), "'"))
  }
  
  if (!is.null(modelBasePath) && modelBasePath != out
      && !dir.create(file.path(modelBasePath, traitFileName), recursive = T)) {
    warning(paste0("Could not create directory '", file.path(modelBasePath, traitFileName), "'"))
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
      message(paste0(
        "Not enough samples present in group '", 
        names(phenotypeFrequencyTable)[phenotypeFrequencyTable < 50], 
        "'. Skipping..."))
      next
    }
    
  } else if (responseDataType == "ordinal") {
    phenotypeFrequencyTable <- table(phenotypeTable$VALUE)
    
    message(paste0("    Available for ", nrow(phenotypeTable), 
                   " samples (classes: ", paste0(names(phenotypeFrequencyTable), collapse = ", "), 
                   ". with the following respective frequencies: ", 
                   paste0(phenotypeFrequencyTable, collapse = ", "), ")."))
    
    if (any(phenotypeFrequencyTable < 50)) {
      message(paste0(
        "Not enough samples present in group '", 
        names(phenotypeFrequencyTable)[phenotypeFrequencyTable < 50], 
        "'. Skipping..."))
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
  
  # Calculate z-score matrix based on polygenic scores and actual phenotypes,
  # using the chosen function for calculating residuals.
  
  pdf(file.path(out, traitFileName, "/debugFigures.pdf"))
  par(xpd = NA)
  
  # Calculate the scaled residuals for every combination
  scaledResidualsMatrix <- calculate.scaledResiduals(
    estimate = completeTable$PGS, 
    actual = completeTable$VALUE, 
    covariates = completeTable[c("AGE", "SEX")],
    responseDataType = responseDataType,
    modelPath = pgsPhenotypeModel)
  
  dev.off()
  
  # The rows correspond to phenotype samples, the cols to genotype samples (polygenic scores)
  rownames(scaledResidualsMatrix) <- completeTable$pheno
  colnames(scaledResidualsMatrix) <- completeTable$geno
  
  if (debug) {
    # Write scaled residuals matrix.
    write.table(scaledResidualsMatrix, file.path(out, traitFileName, "/scaledResidualMatrix.tsv"), 
                sep = "\t", col.names = T, row.names = T, quote = F)
  }

  message("    completed calculating scaled residuals")
  
  message("    calculating log likelihood ratios")
  
  logLikelihoodRatios <- calculate.logLikelihoodRatios(
    valueMatrix = scaledResidualsMatrix, 
    actual = completeTable$VALUE, 
    responseDataType = responseDataType,
    naiveBayesMethod = naiveBayesMethod,
    samplesPerBin = samplesPerNaiveBayesBin,
    classifierPath = modelPath)
  
  rm(scaledResidualsMatrix)
  gc()
  
  if (debug) {
    
    # Write log likelihood ratios.
    write.table(logLikelihoodRatios, file.path(out, traitFileName, "/logLikelihoodRatios.tsv"), 
              sep = "\t", col.names = T, row.names = T, quote = F)
  }
  
  likelihoodRatioDifferenceTest <- t.test(
    diag(logLikelihoodRatios), 
    logLikelihoodRatios[lower.tri(logLikelihoodRatios) | upper.tri(logLikelihoodRatios)],
    alternative = "less")
  
  print(likelihoodRatioDifferenceTest)
  
  if (likelihoodRatioDifferenceTest$p.value <= likelihoodRatioDifferenceAlpha) {
    
    # Resolve log likelihood ratios that are NaN (not a number)
    llrIsNan <- is.nan(logLikelihoodRatios)
    logLikelihoodRatios[llrIsNan] <- 0
    
    # Aggregate number of traits
    aggregatedNumberOfTraits[rownames(logLikelihoodRatios), colnames(logLikelihoodRatios)] <- 
      aggregatedNumberOfTraits[rownames(logLikelihoodRatios), colnames(logLikelihoodRatios)] + !llrIsNan
    
    message("    completed calculating log likelihood ratios.")
    message("    summing log likelihood ratios with aggregated log likelihood ratios...")
    
    # Aggregate log likelihood ratios
    aggregatedLlrMatrix[rownames(logLikelihoodRatios), colnames(logLikelihoodRatios)] <- 
      aggregatedLlrMatrix[rownames(logLikelihoodRatios), colnames(logLikelihoodRatios)] + logLikelihoodRatios
    
    message("    aggregation done!")
  }
  
  llrDataFrame <- 
    as.data.frame.table(logLikelihoodRatios, responseName = "logLikelihoodRatios", stringsAsFactors = FALSE) %>%
    inner_join(link, by = c("Var1" = "pheno"))
  
  rm(logLikelihoodRatios)
  gc()
  
  llrDataFrame$group <- "alternative"
  llrDataFrame$group[llrDataFrame$original == llrDataFrame$Var2] <- "null"
  llrDataFrame$group <- ordered(llrDataFrame$group, levels = c("null", "alternative"))
  
  matrixWideAuc <- auc(
    llrDataFrame$group, llrDataFrame$logLikelihoodRatios)
  
  message(paste0("Calculated overall AUC: ", matrixWideAuc))
  
  permutationTestDataFrame <- llrDataFrame %>%
    filter(geno == Var2)
  
  traitDescriptionsTable[traitIndex, "matrixWideAuc"] <- as.double(matrixWideAuc)
  traitDescriptionsTable[traitIndex, "pValue"] <- likelihoodRatioDifferenceTest$p.value
  
  if (predictingInducedMixUps && "alternative" %in% permutationTestDataFrame$group) {
    
    confinedAuc <- auc(
      permutationTestDataFrame$group, 
      permutationTestDataFrame$logLikelihoodRatios)
    
    message(paste0("Confined AUC: ", confinedAuc))
    traitDescriptionsTable[traitIndex, "confinedAuc"] <- as.double(confinedAuc)
  }

  newAuc <- auc(as.vector(matrix(lower.tri(aggregatedLlrMatrix) | upper.tri(aggregatedLlrMatrix))), as.vector(aggregatedLlrMatrix))
  message(paste0("New AUC: ", newAuc))
  
  # 
  # scaledResidualsDataFrame <-
  #   as.data.frame.table(scaledResidualsMatrix, responseName = "scaledResiduals") %>%
  #   inner_join(link, by = c("Var1" = "pheno")) %>%
  #   mutate(group = case_when(Var2 == geno ~ "null",
  #                            TRUE ~ "alternative"))
  # 
  # scaledResidualsDataFrame$group <- factor(scaledResidualsDataFrame$group, levels = c("null", "alternative"))
  # 
  # rm(scaledResidualsMatrix)
  # gc()
  # 
  # pdf(file.path(out, traitFileName, "unscaledResidualsHistogram.pdf"),
  #     width=8, height = 6, useDingbats = FALSE)
  # 
  # par(xpd = NA)
  # 
  # plotResiduals(scaledResidualsDataFrame, phenotypeTable, responseDataType)
  # 
  # dev.off()
  
  # ggplot(scaledResidualsDataFrame, aes(x=scaledResiduals, stat(density), fill=group)) +
  #   geom_histogram(bins = 72, alpha=.5, position="identity") +
  #   geom_rug(data = scaledResidualsDataFrame %>% filter(geno != original & group == "null"), aes(x=scaledResiduals), inherit.aes=F) +
  #   xlab("Unscaled residuals") + ggtitle(paste0("Unscaled residuals for trait '", trait, "'"))
  
  # ggsave("unscaledResidualsHistogram.png", width=8, height=7)
}

# Write the result values
write.table(traitDescriptionsTable, file.path(out, "outputStatisticsPerTrait.tsv"),
            sep="\t", col.names = T, row.names = T, quote = F)

write.table(aggregatedLlrMatrix, file.path(out, "aggregatedLogLikelihoodRatiosMatrix.tsv"),
            sep="\t", col.names = T, row.names = T, quote = F)

write.table(aggregatedNumberOfTraits, file.path(out, "aggregatedNumberOfTraitsMatrix.tsv"),
            sep="\t", col.names = T, row.names = T, quote = F)

print(dim(aggregatedLlrMatrix))

# aggregatedLlrMatrix <- aggregatedLlrMatrix[apply(aggregatedNumberOfTraits > 0, 1, any), 
#                                            apply(aggregatedNumberOfTraits > 0, 2, any)]

print(dim(aggregatedLlrMatrix))

overallOutputStatistics <- data.frame(name = "overallStatistics")

likelihoodRatioDifferenceTest <- t.test(
  diag(aggregatedLlrMatrix), 
  aggregatedLlrMatrix[lower.tri(aggregatedLlrMatrix) | upper.tri(aggregatedLlrMatrix)],
  alternative = "less")

overallOutputStatistics$pValue <- likelihoodRatioDifferenceTest$p.value

# Convert matrix to data frame for further processing
llrDataFrame <- 
  as.data.frame.table(aggregatedLlrMatrix, responseName = "logLikelihoodRatios") %>%
  inner_join(link, by = c("Var1" = "pheno"))

write.table(llrDataFrame, file.path(out, "/aggregatedLogLikelihoodRatiosDataFrame.tsv"),
            sep="\t", col.names = T, row.names = F, quote = F)

rm(aggregatedLlrMatrix)
gc()

llrDataFrame$group <- "alternative"
llrDataFrame$group[llrDataFrame$original == llrDataFrame$Var2] <- "null"
llrDataFrame$group <- ordered(llrDataFrame$group, levels = c("null", "alternative"))

matrixWideAuc <- auc(
  llrDataFrame$group, llrDataFrame$logLikelihoodRatios)

message(paste0("Calculated overall AUC: ", matrixWideAuc))
overallOutputStatistics$matrixWideAuc <- as.double(matrixWideAuc)

pdf(file.path(out, "ROCcurve_matrixWide.pdf"))
par(xpd = NA)

roc(
  llrDataFrame$group ~ llrDataFrame$logLikelihoodRatios, plot=TRUE,
  print.auc=TRUE,col="green",lwd =4,legacy.axes=TRUE,main="ROC Curves")

dev.off()

# Confine ourselves to the diagonal
permutationTestDataFrame <- llrDataFrame %>%
  filter(geno == Var2)

if ("alternative" %in% permutationTestDataFrame$group) {
  confinedAuc <- auc(
    permutationTestDataFrame$group, 
    permutationTestDataFrame$logLikelihoodRatios)
  
  message(paste0("Confined AUC: ", confinedAuc))
  overallOutputStatistics$confinedAuc <- as.double(confinedAuc)
  
  pdf(file.path(out, "ROCcurve_diagonal.pdf"))
  par(xpd = NA)
  
  roc(
    permutationTestDataFrame$group ~ permutationTestDataFrame$logLikelihoodRatios, 
    plot=TRUE, print.auc=TRUE,col="green",lwd =4,legacy.axes=TRUE,main="ROC Curves")
  
  dev.off()
}

write.table(overallOutputStatistics, file.path(out, "overallOutputStatistics.tsv"),
            sep="\t", col.names = T, row.names = T, quote = F)

ggplot(llrDataFrame, aes(x=logLikelihoodRatios, stat(density), fill=group)) +
  geom_histogram(bins = 32, alpha=.5, position="identity") +
  xlab("Log likelihood ratios") + ggtitle(paste0("LR overall"))

ggsave(file.path(out, "/likelihoodRatioHistogram.png"), width=8, height=7)

ggplot(llrDataFrame, aes(x=group, y=logLikelihoodRatios)) +
  geom_boxplot() + ggtitle(paste0("Log likelihood ratio distributions comparison overall"))

ggsave(file.path(out, "/likelihoodRatioBoxplot.png"), width=8, height=7)

numberOfTraits <- 
  as.data.frame.table(aggregatedNumberOfTraits, responseName = "numberOfTraits")

phenoSamples <- unique(llrDataFrame[llrDataFrame$geno != llrDataFrame$original, "Var1"])

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
