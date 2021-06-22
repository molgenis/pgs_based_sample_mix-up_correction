#!/usr/bin/env Rscript

#############################################################
## c.a.warmerdam@umcg.nl - April 2020
## Processes polygenic scores and phenotypes and
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
                    help='path to a directory to write fitted model parameters to, or to load fitted models from.')
parser$add_argument('--trait-gwas-mapping', required = T,
                    help='path to a tab-delimited file that maps traits to the polygenic scores.')

group <- parser$add_mutually_exclusive_group(required = T)
group$add_argument('--base-pgs-path',
                   help="path to a directory containing polygenic scores. Folder structure should be '<base-pgs-path>/<name-of-gwas-summary-statistic>/<pgs-file-name>'")
group$add_argument('--pgs-file',
                   help="path to a tab-delimited file holding all polygenic score data. Columns should represent sample IDs and individual PGSs respectively.")

parser$add_argument('--pgs-file-name',
                    help="name of files that hold polygenic scores. ignored when a combined table of polygenic scores is supplied (--pgs-file).", 
                    default = "full.UGLI.pgs.profile")

parser$add_argument('--phenotypes-file', required = T,
                    help='path to a tab-delimited file holding all processed phenotype data.')
parser$add_argument('--sample-coupling-file', required = FALSE,
                    help=paste('file containing genotype sample ids in the first column',
                    'and phenotype sample ids in the second column'))
parser$add_argument('--llr-bayes-method', required = FALSE, default = c("NA", "80"), nargs = "+",
                    help = paste('the naive bayes method to use for calculating log likelihood ratios.',
                                  'If NA, the method will be based on the data type for each trait',
                                  '<ewi-discretization | efi-discretization | gaussian | NA> [(average) number of samples per bin]'))

parser$add_argument('--bayes-method-sweep-mode', action='store_true', help = paste0(
  'special mode to perform a sweep over a series of naive bayes methods'))
parser$add_argument('--out',
                    help='path to output directory')
parser$add_argument('--likelihood-ratio-alpha', default = 0.05, help=paste(
  'the t-test alpha to use for selecting traits',
  'based on the difference in log likelihood ratios between the provided and permuted samples'))
parser$add_argument('--output-intermediate-statistics', action='store_true', default = F,
                    help = 'Setting this to false will prevent intermediate AUC calculations, requiring less memory.')

##############################
# Define functions
##############################

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

# Function that returns a matrix of residuals, with phenotype sample across the
# rows, and genotype samples across the columns. If recycle is TRUE, an attempt
# will be made to read this from an existing file:
# - '<intermediateResidualMatrixFileBasePath>/scaledResidualMatrix.tsv', or
# - '<intermediateResidualMatrixFileBasePath>/residualMatrix.rds'.
# If neither of these, exists or can be accessed, the matrix will be created from new.
getResidualsMatrix <- function(
  completeTable, responseDataType, pgsPhenotypeModel, 
  intermediateResidualMatrixFileBasePath, recycle, write) {
  
  intermediateResidualMatrixRdsFilePath <- file.path(
    intermediateResidualMatrixFileBasePath, "residualMatrix.rds")
  intermediateResidualMatrixTsvFilePath <- file.path(
    intermediateResidualMatrixFileBasePath, "scaledResidualMatrix.tsv")
  residualsMatrix <- NULL
  
  if (recycle 
      && file.exists(intermediateResidualMatrixRdsFilePath) 
      && file.access(intermediateResidualMatrixRdsFilePath, 4) == 0) {
    
    message(paste0("    Recycling residuals from '", intermediateResidualMatrixRdsFilePath, "'..."))
    residualsMatrix <- readRDS(intermediateResidualMatrixRdsFilePath)
    
  } else if (recycle 
             && file.exists(intermediateResidualMatrixTsvFilePath) 
             && file.access(intermediateResidualMatrixTsvFilePath, 4) == 0) {
    
    message(paste0("    Recycling residuals from '", intermediateResidualMatrixTsvFilePath, "'..."))
    residualsMatrix <- as.matrix(fread(intermediateResidualMatrixTsvFilePath), rownames = 1)
    
  } else {
    # Calculate residuals matrix based on polygenic scores and actual phenotypes,
    # using the chosen function for calculating residuals.
    
    message("    Calculating residuals...")
    
    # Calculate the scaled residuals for every combination
    residualsMatrix <- calculate.scaledResiduals(
      estimate = completeTable$PGS, 
      actual = completeTable$VALUE, 
      covariates = completeTable[c("AGE", "SEX")],
      responseDataType = responseDataType,
      modelPath = pgsPhenotypeModel)
    
    # The rows correspond to phenotype samples, the cols to genotype samples (polygenic scores)
    rownames(residualsMatrix) <- completeTable$pheno
    colnames(residualsMatrix) <- completeTable$geno
  }
  
  # Check if the residuals matrix is according to what is expected.
  if (!all(dim(residualsMatrix) == c(length(completeTable$pheno), length(completeTable$geno)))) {
    stop("dimensions of residuals matrix do not match the expected dimensions.")
  }
  
  if (!all(rownames(residualsMatrix) == completeTable$pheno)) {
    stop("rownames of residual matrix do not match the expected phenotype sample identifiers.")
  }

  if (!all(colnames(residualsMatrix) == completeTable$geno)) {
    stop("colnames of residual matrix do not match the expected genotype sample identifiers.")
  }
  
  if (write && !file.exists(intermediateResidualMatrixRdsFilePath)) {
    # Write scaled residuals matrix.
    saveRDS(residualsMatrix, intermediateResidualMatrixRdsFilePath)
  }

  return(residualsMatrix)
}

# Function that attempts to creates adapted equal width intervals:
# - 'minFrequencyInTails' values in the outer tails are put in a bin for both tails.
# - the remaining values are split into 'nBins - 2' bins with equal widths.
# - bins that contain less than 'minFrequencyInTails' elements are devided in the centre, 
#   with both sides being added to their closest neighbour bin, until all bins contain at
#   least 'minFrequencyInTails' elements.
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

# Function that takes a series of values, separates them in bins and finds
# likelihood ratios corresponding to each of these bins.
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
gaussianNaiveBayes <- function(values, 
                               parameters) {
  
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
    error = 0.01
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

# Function that returns a matrix of residuals, with phenotype sample across the
# rows, and genotype samples across the columns. If recycle is TRUE, an attempt
# will be made to read this from an existing file:
# - '<intermediateLogLikelihoodRatioMatrixFileBasePath>.logLikelihoodRatios.tsv' (legacy), or
# - '<intermediateLogLikelihoodRatioMatrixFileBasePath>.logLikelihoodRatios.rds'.
# If neither of these, exists or can be accessed, the matrix will be created from new.
getLogLikelihoodRatioMatrix <- function(
  residualsMatrix, completeTable, responseDataType, modelPath, naiveBayesMethod, samplesPerNaiveBayesBin,
  intermediateLogLikelihoodRatioMatrixFileBasePath, recycle, write) {
  
  intermediateLogLikelihoodRatioMatrixRdsFilePath <- paste0(
    intermediateLogLikelihoodRatioMatrixFileBasePath, ".logLikelihoodRatios.rds")
  intermediateLogLikelihoodRatioMatrixTsvFilePath <- paste0(
    intermediateLogLikelihoodRatioMatrixFileBasePath, ".logLikelihoodRatios.tsv")
  LogLikelihoodRatioMatrix <- NULL
  
  if (recycle 
      && file.exists(intermediateLogLikelihoodRatioMatrixRdsFilePath) 
      && file.access(intermediateLogLikelihoodRatioMatrixRdsFilePath, 4) == 0) {
    
    message(paste0("    Recycling log likelihood ratios from '", intermediateLogLikelihoodRatioMatrixRdsFilePath, "'..."))
    LogLikelihoodRatioMatrix <- readRDS(intermediateLogLikelihoodRatioMatrixRdsFilePath)
    
  } else if (recycle 
             && file.exists(intermediateLogLikelihoodRatioMatrixTsvFilePath) 
             && file.access(intermediateLogLikelihoodRatioMatrixTsvFilePath, 4) == 0) {
    
    message(paste0("    Recycling log likelihood ratios from '", intermediateLogLikelihoodRatioMatrixTsvFilePath, "'..."))
    LogLikelihoodRatioMatrix <- as.matrix(fread(intermediateLogLikelihoodRatioMatrixTsvFilePath), rownames = 1)
    
  } else {
    # Calculate residuals matrix based on polygenic scores and actual phenotypes,
    # using the chosen function for calculating residuals.
    
    message("    Calculating residuals...")
    
    # Calculate the scaled residuals for every combination
    LogLikelihoodRatioMatrix <- calculate.logLikelihoodRatios(
      valueMatrix = residualsMatrix,
      actual = completeTable$VALUE,
      responseDataType = responseDataType,
      naiveBayesMethod = naiveBayesMethod,
      samplesPerBin = samplesPerNaiveBayesBin,
      classifierPath = modelPath)
  }
  
  # Check if the residuals matrix is according to what is expected.
  if (!all(dim(LogLikelihoodRatioMatrix) == c(length(completeTable$pheno), length(completeTable$geno)))) {
    stop("dimensions of residuals matrix do not match the expected dimensions.")
  }
  
  if (!all(rownames(LogLikelihoodRatioMatrix) == completeTable$pheno)) {
    stop("rownames of residual matrix do not match the expected phenotype sample identifiers.")
  }
  
  if (!all(colnames(LogLikelihoodRatioMatrix) == completeTable$geno)) {
    stop("colnames of residual matrix do not match the expected genotype sample identifiers.")
  }
  
  if (write && !file.exists(intermediateLogLikelihoodRatioMatrixRdsFilePath)) {
    # Write scaled residuals matrix.
    saveRDS(LogLikelihoodRatioMatrix, intermediateLogLikelihoodRatioMatrixRdsFilePath)
  }
  
  return(LogLikelihoodRatioMatrix)
}

# Function for forcing normal distribution
forceNormal <- function(x) {

  ranked <- rank(x)
  
  pValues <- (0.5 + ranked - 1.0) / length(ranked)
  
  return(qnorm(pValues, 
               mean = mean(x), 
               sd = sd(x)))
}

# Function that adds rows for bayes methods to include in a parameter sweep.
addBayesParameterSweepRows <- function(traitOutputTable) {
  return(traitOutputTable %>%
    select(-naiveBayesMethod, -samplesPerNaiveBayesBin) %>%
    distinct() %>%
    full_join(tibble(naiveBayesMethod = c("gaussian", rep("efi-discretization", 4), rep("ewi-discretization", 4)),
                     samplesPerNaiveBayesBin = c(NA_integer_, rep(c(20, 30, 50, 80), 2))), by = character()) %>%
    mutate(naiveBayesParameters = paste0(naiveBayesMethod, ".", samplesPerNaiveBayesBin)))
}


# Function for loading polygenic scores
loadPolygenicScores <- function(basePathWithPolygenicScores = NULL, pgsFileName = NULL, pgsMergedFile = NULL) {
  
  if (!is.null(basePathWithPolygenicScores)) {
    traitDescriptionsTable$polygenicScoreFilePath <- file.path(basePathWithPolygenicScores, 
                                                               traitDescriptionsTable$summaryStatistics, 
                                                               pgsFileName)
    
    return(bind_rows(mapply(
      function(polygenicScoreFilePath, trait) {
        message(paste0("    Loading polygenic scores from '", polygenicScoreFilePath, "'..."))
        
        # Read the PLINK polygenic score table.
        polygenicScores <- fread(
          polygenicScoreFilePath,
          header=T) %>%
          mutate(TRAIT = trait) %>%
          rename(PGS = SCORESUM, ID = IID) %>%
          select(ID, PGS, TRAIT)
        
        return(polygenicScores)
      },
      traitDescriptionsTable$polygenicScoreFilePath, traitDescriptionsTable$trait, SIMPLIFY = F, USE.NAMES = F)))
    
  } else if (!is.null(pgsMergedFile)) {
    return(fread(pgsMergedFile, header=T, quote="", sep="\t",
                 col.names = c("ID", "TRAIT", "PGS")))
  } else {
    stop("Either 'basePathWithPolygenicScores' or 'pgsMergedFile' must not be NULL.")
  }
}


# Function that makes intermediate plots for residuals
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
args <- parser$parse_args(commandArgs(trailingOnly = TRUE))

message(strwrap(prefix = " ", initial = "", paste(
  "Loading trait-gwas-mapping:\n", args$trait_gwas_mapping)))

# Load table containing paths for the plink output 
# and corresponding phenotype labels.
traitGwasMappingPath <- args$trait_gwas_mapping
traitDescriptionsTable <- fread(
  traitGwasMappingPath, 
  quote="", header=T, sep = "\t", 
  stringsAsFactors=F)

message(strwrap(prefix = " ", initial = "", paste(
  "Loading polygenic scores from:\n", args$trait_gwas_mapping)))

# Get the paths to the polygenic scores.
polygenicScoresTable <- loadPolygenicScores(
  basePathWithPolygenicScores = args$base_pgs_path, pgsFileName = args$pgs_file_name, 
  pgsMergedFile = args$pgs_file)

# Get the output path
out <- args$out

# Should we use debug mode?
debug <- args$debug

# Should we perform a sweep over bayes / likelihood methods for all traits?
loopBayesMethods <- args$bayes_method_sweep_mode

# Output intermediate AUCs?
outputIntermediateStatistics <- args$output_intermediate_statistics

# Set the likelihood alpha
likelihoodRatioDifferenceAlpha <- args$likelihood_ratio_alpha

# Should files that are already present be recycled?
shouldRecycle <- TRUE

# naive bayes / likelihood method
naiveBayesMethod <- args$llr_bayes_method[1]
stopifnot("method must be one of the following: <ewi-discretization | efi-discretization | gaussian>" = 
          is.na(naiveBayesMethod) || naiveBayesMethod %in% c("ewi-discretization", "efi-discretization", "gaussian", "NA"))

# naive bayes / likelihood samples per bin
if (is.na(naiveBayesMethod) 
    || (naiveBayesMethod == "ewi-discretization" 
        || naiveBayesMethod == "efi-discretization" 
        || naiveBayesMethod == "NA") 
    & length(args$llr_bayes_method) == 2) {
  
  samplesPerNaiveBayesBin <- as.numeric(args$llr_bayes_method[2])
}

message(strwrap(prefix = " ", initial = "", paste(
  "Loading phenotype table:\n", args$phenotypes_file)))

# Load the phenotypes 
phenotypesFilePath <- args$phenotypes_file
phenotypesTable <- fread(phenotypesFilePath, header=T, quote="", sep="\t") %>%
  distinct() %>%
  rename_all(recode, "UGLI_ID" = "ID") %>%
  mutate(SEX = case_when(SEX == 1 ~ "Female", SEX == 2 ~ "Male", TRUE ~ as.character(SEX)),
         SEX = factor(SEX, levels = c("Female", "Male"))) %>%
  group_by(ID) %>%
  filter(!any(AGE < 18)) %>%
  ungroup()

link <- data.frame(geno = unique(phenotypesTable$ID), pheno = unique(phenotypesTable$ID), stringsAsFactors = F) %>%
  slice_sample(10000)

# Get the link path
if (!is.null(args$sample_coupling_file)) {
  message(strwrap(prefix = " ", initial = "", paste(
    "Loading sample coupling table:\n", args$sample_coupling_file)))
  
  sampleCouplingFilePath <- args$sample_coupling_file
  link <- fread(sampleCouplingFilePath, stringsAsFactors=F, header=T) %>%
    filter(pheno %in% unique(phenotypesTable$ID))
} else {
  message(strwrap(prefix = " ", initial = "", 
                  "Generating sample coupling table from phenotypes file."))
}

predictingInducedMixUps <- F

if (!("original" %in% colnames(link))) {
  link$original <- link$geno
} else {
  predictingInducedMixUps <- T
}

# Link the phenotypes table to the sample coupling file
phenotypesTable <- phenotypesTable %>%
  rename(pheno = ID) %>%
  inner_join(link, by="pheno")

# Merge the PGSs with actual phenotype data
completeLongTable <- phenotypesTable %>%
  inner_join(polygenicScoresTable, by = c("geno" = "ID", "TRAIT" = "TRAIT"))

modelBasePath <- NULL

if (!is.null(args$base_fit_model_path)) {
  modelBasePath <- args$base_fit_model_path
}

aggregatedLlrMatrix <- matrix(nrow = nrow(link), 
                              ncol = nrow(link), 
                              dimnames = list(link$pheno, link$geno), 0)

aggregatedNumberOfTraits <- matrix(nrow = nrow(link), 
                                   ncol = nrow(link), 
                                   dimnames = list(link$pheno, link$geno), 0)

traitDescriptionsTable <- traitDescriptionsTable %>%
  group_by(trait) %>%
  mutate(naiveBayesMethod = case_when(
    !is.na(naiveBayesMethod) 
    & naiveBayesMethod %in% c("gaussian", "efi-discretization", "ewi-discretization") ~ as.character(naiveBayesMethod),
    traitDataType == "continuous" ~ "gaussian",
    traitDataType %in% c("ordinal", "binary") ~ "ewi-discretization"),
    samplesPerNaiveBayesBin = samplesPerNaiveBayesBin,
    naiveBayesParameters = paste0(naiveBayesMethod, ".", samplesPerNaiveBayesBin))

traitOutputTable <- traitDescriptionsTable %>%
  mutate(confinedAuc = NA_real_,
         confinedAucOnScaledLlr = NA_real_,
         matrixWideAuc = NA_real_,
         pValue = NA_real_,
         traitOutputDir = NA_character_,
         modelBasePath = modelBasePath)

if (loopBayesMethods) {
  traitOutputTable <- addBayesParameterSweepRows(traitOutputTable)
}

# Loop trough traits

for (traitIndex in 1:nrow(traitDescriptionsTable)) {
  
  trait <- traitDescriptionsTable$trait[traitIndex]
  responseDataType <- traitDescriptionsTable$traitDataType[traitIndex]
  
  traitFileName <- paste(traitIndex, gsub(" ", "_", trait), sep = ".")
  traitDirectory <- file.path(out, traitFileName)
  traitOutputTable[traitOutputTable$trait == trait, "traitOutputDir"] <- traitFileName
  
  modelPath <- NULL
  pgsPhenotypeModel <- NULL
  
  if (!is.null(modelBasePath)) {
    pgsPhenotypeModel <- file.path(modelBasePath, traitFileName, "pgsPhenotypeModel.rda")
  }

  if (!dir.create(traitDirectory, recursive = T, showWarnings = FALSE)) {
    warning(paste0("Could not create directory '", traitDirectory, "'"))
  }
  
  if (!is.null(modelBasePath) && modelBasePath != out
      && !dir.create(file.path(modelBasePath, traitFileName), recursive = T, showWarnings = FALSE)) {
    warning(paste0("Could not create directory '", file.path(modelBasePath, traitFileName), "'"))
  }
  
  # Filter the completeLongTable table
  completeTable <- completeLongTable %>%
    filter(TRAIT == trait)
  
  traitOutputTable[traitOutputTable$trait == trait, "numberOfSamples"] <- nrow(completeTable)
  
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
    message(paste0("    Available for ", nrow(completeTable), " samples."))
  }

  # Calculate residuals matrix
  residualsMatrix <- getResidualsMatrix(
    completeTable, responseDataType, pgsPhenotypeModel, 
    traitDirectory, recycle = shouldRecycle, write = debug)

  message("    completed calculating scaled residuals")
  
  for (naiveBayesParameters in traitOutputTable$naiveBayesParameters[traitOutputTable$trait == trait]) {
    
    naiveBayesMethod <- traitOutputTable[traitOutputTable$trait == trait & traitOutputTable$naiveBayesParameters == naiveBayesParameters, 
                                         "naiveBayesMethod"]
    samplesPerNaiveBayesBin <- traitOutputTable[traitOutputTable$trait == trait & traitOutputTable$naiveBayesParameters == naiveBayesParameters, 
                                         "samplesPerNaiveBayesBin"]
    
    message(paste0("    log likelihood method: ", naiveBayesMethod, ", ", samplesPerNaiveBayesBin))
    
    if (!is.null(modelBasePath)) {
      modelPath <- file.path(modelBasePath, traitFileName, naiveBayesParameters)
      if (!dir.create(modelPath, recursive = T)) {
        warning(paste0("Could not create directory '", modelPath, "'"))
      }
    }
    
    intermediateLogLikelihoodRatioMatrixFileBasePath <- file.path(
      traitDirectory, naiveBayesParameters)
    
    logLikelihoodRatios <- getLogLikelihoodRatioMatrix(
      residualsMatrix, completeTable, responseDataType, 
      modelPath, naiveBayesMethod, samplesPerNaiveBayesBin,
      intermediateLogLikelihoodRatioMatrixFileBasePath, recycle = shouldRecycle, write  = debug)
    
    likelihoodRatioDifferenceTest <- t.test(
      diag(logLikelihoodRatios), 
      logLikelihoodRatios[lower.tri(logLikelihoodRatios) | upper.tri(logLikelihoodRatios)],
      alternative = "less")
    
    print(likelihoodRatioDifferenceTest)
    
    if (likelihoodRatioDifferenceTest$p.value <= likelihoodRatioDifferenceAlpha & !loopBayesMethods) {
      
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
    
    traitOutputTable[traitOutputTable$trait == trait & traitOutputTable$naiveBayesParameters == naiveBayesParameters, "pValue"] <- 
      likelihoodRatioDifferenceTest$p.value
    
    # If output of intermediate statistics is requested, calculate these.
    if (outputIntermediateStatistics) {
    
      llrDataFrame <- 
        as.data.frame.table(logLikelihoodRatios, responseName = "logLikelihoodRatios", stringsAsFactors = FALSE) %>%
        inner_join(link, by = c("Var1" = "pheno"))
      
      rm(logLikelihoodRatios)
      gc()
      
      llrDataFrame$group <- "alternative"
      llrDataFrame$group[llrDataFrame$original == llrDataFrame$Var2] <- "null"
      llrDataFrame$group <- ordered(llrDataFrame$group, levels = c("null", "alternative"))
      
      llrDataFrame <- llrDataFrame %>% 
        group_by(Var1) %>%
        mutate(scaledLlr = scale(logLikelihoodRatios)[,1])
      
      matrixWideAucOnScaledLlr <- auc(
        llrDataFrame$group, llrDataFrame$scaledLlr)
      
      message(paste0("Calculated overall AUC on scaled log likelihood ratios: ", matrixWideAucOnScaledLlr))
      
      permutationTestDataFrame <- llrDataFrame %>%
        filter(geno == Var2)
      
      rm(llrDataFrame)
      gc()
      
      traitOutputTable[traitOutputTable$trait == trait & traitOutputTable$naiveBayesParameters == naiveBayesParameters, "matrixWideAucOnScaledLlr"] <- 
        as.double(matrixWideAucOnScaledLlr)
      
      if (predictingInducedMixUps && "alternative" %in% permutationTestDataFrame$group) {
        
        confinedAuc <- auc(
          permutationTestDataFrame$group, 
          permutationTestDataFrame$logLikelihoodRatios)
        
        confinedAucOnScaledLlr <- auc(
          permutationTestDataFrame$group, 
          permutationTestDataFrame$scaledLlr)
        
        message(paste0("Confined AUC: ", confinedAuc))
        traitOutputTable[traitOutputTable$trait == trait & traitOutputTable$naiveBayesParameters == naiveBayesParameters, "confinedAuc"] <- 
          as.double(confinedAuc)
        
        message(paste0("Confined AUC on scaled log likelihood ratios: ", confinedAucOnScaledLlr))
        traitOutputTable[traitOutputTable$trait == trait & traitOutputTable$naiveBayesParameters == naiveBayesParameters, "confinedAucOnScaledLlr"] <- 
          as.double(confinedAucOnScaledLlr)
      }
      
      rm(permutationTestDataFrame)
      gc()
    }
    
    # Clear the intermediate residual and log likelihood ratio matrix
    rm(residualsMatrix)
    rm(logLikelihoodRatios)
    gc()
  }
}

# Write the result values
message(paste0("Exporting output statistics per trait: 'outputStatisticsPerTrait.tsv'"))
write.table(traitOutputTable, file.path(out, "outputStatisticsPerTrait.tsv"),
            sep="\t", col.names = T, row.names = T, quote = F)

message(paste0("Exporting raw aggregated log likelihood ratio matrix: 'aggregatedLogLikelihoodRatiosMatrix.tsv'"))
write.table(aggregatedLlrMatrix, file.path(out, "aggregatedLogLikelihoodRatiosMatrix.tsv"),
            sep="\t", col.names = T, row.names = T, quote = F)

message(paste0("Exporting number of traits for every phenotype-genotype combination: 'aggregatedNumberOfTraitsMatrix.tsv'"))
write.table(aggregatedNumberOfTraits, file.path(out, "aggregatedNumberOfTraitsMatrix.tsv"),
            sep="\t", col.names = T, row.names = T, quote = F)

scaledLogLikelihoodMatrix <- t(apply(aggregatedLlrMatrix, 1, function(x) scale(x)))

# Extract the diagonal values
diagValues <- aggregatedLlrMatrix[lower.tri(aggregatedLlrMatrix, diag = TRUE)
                                  & upper.tri(aggregatedLlrMatrix, diag = TRUE)]

# Extract the diagonal values
diagValuesScaled <- scaledLogLikelihoodMatrix[lower.tri(scaledLogLikelihoodMatrix, diag = TRUE)
                                              & upper.tri(scaledLogLikelihoodMatrix, diag = TRUE)]

# Extract the diagonal values
diagTraitNumbers <- aggregatedNumberOfTraits[lower.tri(aggregatedNumberOfTraits, diag = TRUE)
                                             & upper.tri(aggregatedNumberOfTraits, diag = TRUE)]

results <- tibble(Var1 = rownames(aggregatedLlrMatrix),
                  Var2 = colnames(aggregatedLlrMatrix),
                  logLikelihoodRatios = diagValues,
                  scaledLlr = diagValuesScaled, 
                  numberOfTraits = diagTraitNumbers)

message(paste0("Exporting output matrix: 'idefixPredictions.txt'"))

# Write Idéfix predictions.
write.table(results, file.path(out, "idefixPredictions.txt"), row.names = F, col.names = T, quote = F, sep = "\t")

if (loopBayesMethods) {
  stop("exiting...")
}

# Define a table to write statistics to
overallOutputStatistics <- data.frame(name = "overallStatistics")

# Convert matrix to data frame for further processing
llrDataFrame <- 
  as.data.frame.table(aggregatedLlrMatrix, responseName = "logLikelihoodRatios") %>%
  inner_join(link, by = c("Var1" = "pheno")) %>%
  mutate(
    diag = case_when(
      geno == Var2 ~ T,
      geno != Var2 ~ F),
    correct = case_when(
      original == Var2 ~ T,
      original != Var2 ~ F),
    mixUp = case_when(
      diag & correct ~ F,
      diag & !correct ~ T),
    group = case_when(diag & !correct ~ "inducedMixUp",
                      diag & correct ~ "provided",
                      !diag ~ "permuted")) %>%
  group_by(Var1) %>%
  mutate(scaledLlr = scale(logLikelihoodRatios)[,1])

# Remove the aggregated llr matrix in favour of the data frame
rm(aggregatedLlrMatrix)
gc()

numberOfTriatsDiag <- diag(aggregatedNumberOfTraits)

rm(aggregatedNumberOfTraits)
gc()

matrixWideAucOnScaledLlr <- auc(
  llrDataFrame$correct, 
  llrDataFrame$scaledLlr)

message(paste0("Matrix-wide AUC on scaled log likelihood ratios: ", matrixWideAucOnScaledLlr))
overallOutputStatistics$matrixWideAuc <- as.double(matrixWideAucOnScaledLlr)

message(paste0("Exporting matrix-wide ROC curve: 'ROCcurve_matrixWide_scaled.pdf'"))

pdf(file.path(out, "ROCcurve_matrixWide_scaled.pdf"))
par(xpd = NA)

# Calculate plot
roc(
  llrDataFrame$mixUp ~ llrDataFrame$scaledLlr, plot=TRUE, 
  print.auc=TRUE,col="green",lwd =4,legacy.axes=TRUE,main="ROC Curves")

dev.off()

# Confine ourselves to the diagonal
permutationTestDataFrame <- llrDataFrame %>%
  filter(geno == Var2)

rm(llrDataFrame)
gc()

permutationTestDataFrame$numberOfTraits <- numberOfTriatsDiag

message(paste0("Exporting output matrix: 'providedSampleDataFrame.tsv'"))

# Export the log likelihood data frame with scaled values
write.table(permutationTestDataFrame, file.path(out, "providedSampleDataFrame.tsv"),
            sep="\t", col.names = T, row.names = F, quote = F)

# Calculate the mix-ups
if ("inducedMixUp" %in% permutationTestDataFrame$group) {
  
  confinedAucOnScaledLlr <- auc(
    permutationTestDataFrame$group, 
    permutationTestDataFrame$scaledLlr)
    
  message(paste0("Confined AUC on scaled log likelihood ratios: ", confinedAucOnScaledLlr))
  overallOutputStatistics$confinedAuc <- as.double(confinedAucOnScaledLlr)
  
  message(paste0("Exporting ROC curve for provided samples: 'ROCcurve_diagonal_scaled.pdf'"))
  
  pdf(file.path(out, "ROCcurve_diagonal_scaled.pdf"))
  par(xpd = NA)

  roc(
    permutationTestDataFrame$group ~ permutationTestDataFrame$scaledLlr, 
    plot=TRUE, print.auc=TRUE,col="green",lwd =4,legacy.axes=TRUE,main="ROC Curves on scaled LLR")
  
  dev.off()
  
  ggplot(permutationTestDataFrame, aes(x=scaledLlr, stat(density), fill=group)) +
    geom_histogram(bins = 32, alpha=.5, position="identity") +
    xlab("Log likelihood ratios") + ggtitle(paste0("scaled LLR overall"))
  
  ggsave(file.path(out, "scaledlikelihoodRatioHistogram.png"), width=8, height=7)
}

write.table(overallOutputStatistics, file.path(out, "overallOutputStatistics.tsv"),
            sep="\t", col.names = T, row.names = T, quote = F)

message("Done!")
