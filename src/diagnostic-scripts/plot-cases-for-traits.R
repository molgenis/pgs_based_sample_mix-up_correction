#!/usr/bin/env Rscript

#############################################################
## c.a.warmerdam@umcg.nl - November 2020
## Processes plink polygenic scores and phenotypes and
## calculates predictive power for every trait
#############################################################

##############################
# Load libraries
##############################
library(MASS)
library(pROC)
library(tidyverse)
library(argparse)
library(data.table)
library(stringr)
library(ggrepel)
library(sfsmisc)

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
parser$add_argument('--out',
                    help='path to output directory')

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
  
  return(residuals)
}

# Function that returns a function for calculating residuals.
# The model determining residuals is dependent on the type of response data type.
responseFunConstructor <- function(estimate, actual, covariates, 
                                   modelPath, useOrderedLogit = TRUE) {
  
  # If the path to write/read fitted models to already points to an existing file, load this.
  # In this case residual calculation should be performed using the loaded model
  load(modelPath)
  
  # Return a linear model in case 'actual', is a continuous data type.
  if (exists("olsModel")) {

    print(summary(olsModel))
    
    # Define residuals function based on the linear regression model.
    responseFun <- function(estimate, actual, covariates) {
      return(responseFromOlsRegressionLine(estimate = estimate, 
                                           actual = actual, 
                                           covariates = covariates, 
                                           olsModel = olsModel))
    }
    
    return(responseFun)
    # Return a logistic model in case 'actual', is a binary data type.
  } else if (exists("logitModel")) {
    
    print(summary(logitModel))
    
    # Define residuals function using the logistic regression model.
    responseFun <- function(estimate, actual, covariates) {
      return(responseFromLogitRegressionLine(estimate = estimate, 
                                             actual = actual, 
                                             covariates = covariates, 
                                             logitModel = logitModel))
    }
    
    return(responseFun)
    
    # Return residuals function based on ordered logit in case 'actual' is an ordinal data type,
    # other than binary
  } else if (exists("orderedLogitModel")) {
    
    classes <- sort(unique(actual))
    actualOrdered <- factor(actual, levels = classes, ordered = TRUE)
    
    print(summary(orderedLogitModel))
    
    # Define residuals function using the ordered logistic regression model.
    responseFun <- function(estimate, actual, covariates) {
      actualOrdered <- factor(actual, levels = classes, ordered = TRUE)
      
      responseProbabilities <- responseFromOrderedLogitModel(
        estimate = estimate, actual = actualOrdered,
        covariates = covariates, orderedLogitModel = orderedLogitModel)
      
      classes * responseProbabilities
    }
    
    return(responseFun)
  }
}

calculateResponseValues <- function(estimate, actual, covariates, modelPath, useOrderedLogit = TRUE) {
  responseFun <- responseFunConstructor(estimate = estimate,
                                        actual = actual,
                                        covariates = covariates,
                                        modelPath = modelPath,
                                        useOrderedLogit = useOrderedLogit)
  
  # extract residuals.
  response <- responseFun(estimate = estimate, 
                            actual = actual, 
                            covariates = covariates)
  
  return(response)
}

responseFromOlsRegressionLine <- function(estimate, actual, covariates, olsModel) {
  # Predict actual values using the given model and the supplied independent variables
  predictedActual <- predict(olsModel,
                             cbind(estimate, covariates), 
                             type = "response")
  
  return(unname(predictedActual))
}

responseFromLogitRegressionLine <- function(estimate, actual, covariates, logitModel) {
  # Predict actual values using the given model and the supplied independent variables.
  predictedActual <- predict(logitModel, 
                             cbind(estimate, covariates), 
                             type = "response")
  
  return(unname(predictedActual))
}

responseFromOrderedLogitModel <- function(estimate, actual, covariates, orderedLogitModel) {
  # Predict actual values using the given model and the supplied independent variables.
  predictedProbabilities <- predict(orderedLogitModel, 
                                    cbind(estimate, covariates), 
                                    type = "probs")

  return(predictedProbabilities[cbind(
    1:nrow(predictedProbabilities), 
    actual)])
}

# Define function for calculating the AUC
calculate.auc <- function(actual, predictor) {
  pr <- prediction(predictor, actual)
  
  performance <- performance(pr, measure = "auc")
  auc <- performance@y.values[[1]]
  return(auc)
}

read.diagonal = function(filepath, delim = "\t") {
  
  colNames <- c()
  rowNames <- c()
  values <- c()
  
  # Open the file
  con <- file(filepath, "r")
  
  # Start at the first line
  index <- 0
  
  # The last line should have no length, break if this happens
  while ( TRUE ) {
    index <- index + 1
    
    # Read the first line
    line <- readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    
    # Process individual lines
    # The first line should hold column names
    if (index == 1) {
      splitted_string <- str_split(line, fixed(delim))[[1]]
      colNames <- splitted_string
    }
    # The other lines should hold the value on interest on the correct column
    else if (index != 1) {
      splitted_string <- str_split(line, fixed(delim), n = index + 1)[[1]]
      rowNames <- c(rowNames, splitted_string[1])
      values <- c(values, as.double(splitted_string[index]))
    }
  }
  
  close(con)
  
  stopifnot(!is.null(colNames) && !is.null(rowNames) && colNames == rowNames)
  names(values) <- colNames
  return(values)
}

read.samplesOfInterest = function(filepath, samples, delim = "\t") {
  
  colNames <- c()
  rowNames <- c()
  values <- c()
  
  # Open the file
  con <- file(filepath, "r")
  
  # Start at the first line
  index <- 0
  
  # The last line should have no length, break if this happens
  while ( TRUE ) {
    index <- index + 1
    
    # Read the first line
    line <- readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    
    # Process individual lines
    # The first line should hold column names
    if (index == 1) {
      splitted_string <- str_split(line, fixed(delim))[[1]]
      colNames <- splitted_string
    }
    # The other lines should hold the value on interest on the correct column
    else if (index != 1 && colNames[index - 1] %in% samples) {
      splitted_string <- str_split(line, fixed(delim))[[1]]
      rowNames <- c(rowNames, splitted_string[1])
      values <- c(values, as.double(splitted_string[2:length(splitted_string)]))
    }
    cat('\r',as.character(index))
    flush.console() 
  }
  
  close(con)
  
  dim(values) <- c(length(colNames), length(rowNames))
  values <- t(values)
  dimnames(values) <- list(rowNames, colNames)
  
  return(values)
}

read.matrix = function(filepath, delim = "\t") {
  
  colNames <- c()
  rowNames <- c()
  values <- c()
  
  # Open the file
  con <- file(filepath, "r")
  
  # Start at the first line
  index <- 0
  
  # The last line should have no length, break if this happens
  while ( TRUE ) {
    index <- index + 1
    
    # Read the first line
    line <- readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    
    # Process individual lines
    # The first line should hold column names
    if (index == 1) {
      splitted_string <- str_split(line, fixed(delim))[[1]]
      colNames <- splitted_string
    }
    # The other lines should hold the value on interest on the correct column
    else if (index != 1) {
      splitted_string <- str_split(line, fixed(delim))[[1]]
      rowNames <- c(rowNames, splitted_string[1])
      values <- c(values, as.double(splitted_string[2:length(splitted_string)]))
    }
    cat('\r',as.character(index))
    flush.console() 
  }
  
  close(con)
  
  dim(values) <- c(length(colNames), length(rowNames))
  values <- t(values)
  dimnames(values) <- list(rowNames, colNames)
  
  return(values)
}

getLogLikelihoodRatiosForSoi <- function(trait, traitIndex, samples) {
  
  polygenicScoreFilePath <- traitDescriptionsTable$polygenicScoreFilePath[traitIndex]
  stopifnot(traitDescriptionsTable$trait[traitIndex] == trait)
  responseDataType <- traitDescriptionsTable$traitDataType[traitIndex]
  
  traitFileName <- paste(traitIndex, gsub(" ", "_", trait), sep = ".")
  
  naiveBayesParameters <- traitDescriptionsTable[traitDescriptionsTable$trait == trait, "naiveBayesParameters"]
  
  intermediateLogLikelihoodRatiosFilePath <- file.path(
    out, traitFileName, paste0(naiveBayesParameters, ".", "logLikelihoodRatios.tsv"))
  
  if (!(file.exists(intermediateLogLikelihoodRatiosFilePath) 
        && file.access(intermediateLogLikelihoodRatiosFilePath, 4) == 0)) {
    return()
  }
  
  message(paste0("    Recycling log likelihood ratios from '", intermediateLogLikelihoodRatiosFilePath, "'..."))
  # table_tmp <- read_delim(intermediateLogLikelihoodRatiosFilePath, delim = "\t", quote = "")
  # logLikelihoodRatios <- as.matrix(table_tmp[,-1])
  # rownames(logLikelihoodRatios) <- table_tmp[,1][[1]]
  # rm(table_tmp)
  
  
  logLikelihoodRatiosDataFrame <- logLikelihoodRatios %>%
    as.data.frame.table(responseName = "logLikelihoodRatios", stringsAsFactors = FALSE) %>%
    group_by(Var1) %>%
    mutate(scaledLogLikelihoodRatios = scale(logLikelihoodRatios)[,1]) %>%
    ungroup() %>%
    filter(Var2 %in% samples)
  
  rm(logLikelihoodRatios)
  gc()
  
  return(logLikelihoodRatiosDataFrame)
}

scalePerGroup <- function(values, group) {
  return(tibble(values = values, grouping = group) %>%
    group_by(grouping) %>%
    mutate(scaled = scale(values)[,1]) %>%
    ungroup() %>%
    pull(scaled))
}

forceNormal <- function(x) {
  
  ranked <- rank(x)
  
  pValues <- (0.5 + ranked - 1.0) / length(ranked)
  
  return(qnorm(pValues, 
               mean = mean(x), 
               sd = sd(x)))
}


##############################
# Run
##############################
args <- parser$parse_args(c("--trait-gwas-mapping", "/groups/umcg-lld/tmp01/other-users/umcg-rwarmerdam/pgs_based_mixup_correction/scripts/r-scripts/pgs_based_sample_mix-up_correction/data/lifelines/trait-gwas-mapping.txt",
                            "--base-pgs-path", "/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/output/PRScs/20201120/",
                            "--phenotypes-file", "/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/data/lifelines/processed/pgs.phenotypes_20201215.ugli.dat",
                            "--out", "/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/output/sample-swap-prediction/20201120.20210127.NA.80/"))

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
basePathWithPolygenicScores <- args$base_pgs_path
traitDescriptionsTable$polygenicScoreFilePath <- file.path(basePathWithPolygenicScores, 
                                                           traitDescriptionsTable$summaryStatistics, 
                                                           "full.UGLI.pgs.profile")

# Get the output path
out <- args$out
modelBasePath <- file.path(out, "case-visualization-models")

if (!dir.create(out, recursive = T)) {
  warning(paste0("Could not create directory '", out, "'"))
}

message(strwrap(prefix = " ", initial = "", paste(
  "Loading phenotype table:\n", args$phenotypes_file)))

# Load the phenotypes 
phenotypesFilePath <- args$phenotypes_file
phenotypesTable <- fread(phenotypesFilePath, header=T, quote="", sep="\t") %>%
  rename_all(recode, "UGLI_ID" = "ID") %>%
  mutate(SEX = factor(SEX, levels = c("Female", "Male"))) %>%
  group_by(ID) %>%
  filter(!any(AGE < 18) & !is.na(VALUE)) %>%
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

qc <- fread("/groups/umcg-lifelines/tmp01/releases/gsa_genotypes/v1/Logs/ugli_qc_release1_samples.csv")
gsaLink <- fread("/groups/umcg-lifelines/tmp01/releases/gsa_linkage_files/v1/gsa_linkage_file.dat")
idefixPredictions <- fread(file.path("/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/output/sample-swap-prediction/20201120.20210127.NA.80/", "UGLI_idefix_predictions.txt")) %>%
  inner_join(gsaLink, by = c("UGLI_ID")) %>%
  inner_join(qc, by = c("genotyping_name" = "Sample_ID"))

nSamplesToSelect <- 6
nSamplesToSelectGood <- 4
nSamplesToSelectBad <- 4

sampleTable <- bind_rows(idefixPredictions %>% slice_max(scaledLogLikelihoodRatios, n = nSamplesToSelectBad) %>%
                           arrange(desc(scaledLogLikelihoodRatios)) %>%
                           mutate(sampleGoodSampleBad = factor(c("bad"), levels = c("good", "neutral", "bad")),
                                  index = 1:nSamplesToSelectBad),
                         idefixPredictions %>% slice_min(scaledLogLikelihoodRatios, n = nSamplesToSelectGood) %>%
                           arrange(scaledLogLikelihoodRatios) %>%
                           mutate(sampleGoodSampleBad = factor(c("good"), levels = c("good", "neutral", "bad")),
                                  index = 1:nSamplesToSelectGood)) %>%
  rename(ID = UGLI_ID) %>%
  dplyr::select(ID, sampleGoodSampleBad, index, scaledLogLikelihoodRatios)

traitDescriptionsTable <- traitDescriptionsTable %>%
  group_by(trait) %>%
  mutate(naiveBayesMethod = case_when(
    traitDataType == "continuous" ~ "gaussian",
    traitDataType %in% c("ordinal", "binary") ~ "ewi-discretization"),
    samplesPerNaiveBayesBin = 80,
    naiveBayesParameters = paste0(naiveBayesMethod, ".", samplesPerNaiveBayesBin))

logLikelihoodRatiosPerTrait <- mapply(getLogLikelihoodRatiosForSoi, 
                                      traitDescriptionsTable$trait, 
                                      1:nrow(traitDescriptionsTable),
                                      MoreArgs = list(samples = sampleTable$ID),
                                      SIMPLIFY = F, USE.NAMES = T)

logLikelihoodRatiosPerTraitCombined <- bind_rows(logLikelihoodRatiosPerTrait, .id = "TRAIT") %>%
  inner_join(sampleTable, by = c("Var1" = "ID")) %>%
  filter(Var1 == Var2)

#saveRDS(logLikelihoodRatiosPerTraitCombined, file.path(out, "logLikelihoodRatiosPerTrait.rds"))

pdf(file.path(out, paste0("logLikelihoodRatiosPerSoi_", format(Sys.Date(), "%Y%m%d"), ".pdf")), 
    useDingbats = FALSE, width = 8, height = 11.5)

par(xpd = NA)

ggplot(data = logLikelihoodRatiosPerTraitCombined,
       aes(y = scaledLogLikelihoodRatios.x, x = TRAIT, colour = sampleGoodSampleBad)) +
  geom_col() +
  coord_flip() +
  facet_grid(index ~ sampleGoodSampleBad)

# If you map an aesthetic to a categorical variable, you will get a
# set of contours for each value of that variable

dev.off()

completeTables <- list()


  
# Loop trough traits
for (traitIndex in 1:nrow(traitDescriptionsTable)) {
  
  polygenicScoreFilePath <- traitDescriptionsTable$polygenicScoreFilePath[traitIndex]
  trait <- traitDescriptionsTable$trait[traitIndex]
  responseDataType <- traitDescriptionsTable$traitDataType[traitIndex]
  
  traitFileName <- paste(traitIndex, gsub(" ", "_", trait), sep = ".")

  pgsPhenotypeModel <- NULL
  
  if (!is.null(modelBasePath)) {
    pgsPhenotypeModel <- file.path(modelBasePath, traitFileName, "pgsPhenotypeModel.rda")
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
    
  } else if (responseDataType == "ordinal") {
    phenotypeFrequencyTable <- table(phenotypeTable$VALUE)
    
    message(paste0("    Available for ", nrow(phenotypeTable), 
                   " samples (classes: ", paste0(names(phenotypeFrequencyTable), collapse = ", "), 
                   ". with the following respective frequencies: ", 
                   paste0(phenotypeFrequencyTable, collapse = ", "), ")."))
    
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
  
  residuals <- calculate.scaledResiduals(
    estimate = completeTable$PGS, 
    actual = completeTable$VALUE, 
    covariates = completeTable[c("AGE", "SEX")],
    responseDataType = responseDataType,
    modelPath = pgsPhenotypeModel)
  
  predictedPhenotypes <- calculateResponseValues(
    estimate = completeTable$PGS, 
    actual = completeTable$VALUE, 
    covariates = completeTable[c("AGE", "SEX")],
    modelPath = pgsPhenotypeModel)
  
  # Get the residual matrix:
  # Either from an existing output directory, or
  # Calculate it again here.
  
  # We only need the diagonal
  # residuals <- read.diagonal(file.path(out, traitFileName, "scaledResidualMatrix.tsv"))
  #residuals <- read.diagonal("test123.txt")

  completeTables[[trait]] <- tibble(ID = completeTable$pheno, 
                                    residuals = residuals,
                                    predictions = predictedPhenotypes,
                                    ) %>%
    inner_join(completeTable, by = c("ID" = "pheno"))
}

completeTablesCombined <- bind_rows(completeTables, .id = "TRAIT") %>%
  ungroup() %>%
  dplyr::select(ID, PGS, VALUE, TRAIT, 
                residuals, predictions) %>%
  left_join(sampleTable, by = c("ID")) %>%
  left_join(traitDescriptionsTable %>% dplyr::select(trait, traitDataType), by = c("TRAIT" = "trait")) %>%
  mutate(sampleGoodSampleBad = case_when(
    !is.na(sampleGoodSampleBad) ~ sampleGoodSampleBad,
    TRUE ~ factor("neutral", levels = c("good", "neutral", "bad")))) %>%
  group_by(TRAIT) %>%
  mutate(
    actual = predictions + residuals,
    scaledResiduals = case_when(
      traitDataType != "continuous" ~ scalePerGroup(residuals, VALUE),
      TRUE ~ scale(residuals)[,1]),
    scaledActual = case_when(
      traitDataType != "continuous" ~ scalePerGroup(actual, VALUE),
      TRUE ~ scale(actual)[,1]),
    scaledPredictions = forceNormal(scale(predictions)[,1]),
    scaledAltActual = scaledPredictions + scaledResiduals
  ) %>%
  ungroup() %>%
  arrange(scaledLogLikelihoodRatios) %>%
  mutate(pLog10 = -log10(2*pnorm(-abs(scaledResiduals))),
         traitLabel = case_when(abs(scaledResiduals) > 2.58 ~ paste0(gsub(" ", "~", TRAIT), "~(italic(p) == ", pretty10exp(2*pnorm(-abs(scaledResiduals)), digits = 2), ")"),
                                TRUE ~ ""),
         traitLabelSmall = case_when(abs(scaledResiduals) > 2.58 ~ paste0("italic(p) == ", pretty10exp(2*pnorm(-abs(scaledResiduals)), digits = 2)),
                                TRUE ~ ""),
         sampleLabel = paste0(str_to_title(sampleGoodSampleBad), "~sample~italic(", index, ")~(prediction == ", round(scaledLogLikelihoodRatios, digits = 1), ")"),
         sampleLabel = ordered(sampleLabel, levels = unique(sampleLabel)),
         indexLabel = paste0(index, ":~(prediction == ", round(scaledLogLikelihoodRatios, digits = 1), ")"),
         sampleGoodSampleBadLabel = case_when(sampleGoodSampleBad == "good" ~ "Correct sample",
                                              sampleGoodSampleBad == "bad" ~ "Predicted mix-up")) %>%
  group_by(ID) %>%
  mutate(overallZscore = sum(scaledResiduals) / sqrt(n())) %>%
  ungroup()

# completeTablesSummarised <- completeTablesCombined %>%
#   filter(traitDataType == "continuous") %>%
#   group_by(TRAIT) %>%
#   summarise(mean = mean(scaledResiduals),
#             sd = sd(scaledResiduals))

completeTablesSummarised <- completeTablesCombined %>%
  group_by(ID) %>%
  mutate(overallZscore = sum(scaledResiduals) / sqrt(n()))

old <- theme_set(theme_classic())
theme_update(line = element_line(
  colour = "black", size = (0.5 * 72.27/96), 
  linetype = 1, lineend = "butt", arrow = F, inherit.blank = T),
  strip.background = element_rect(colour = NA, fill = NA),
  axis.line = element_line(
    colour = "#595A5C", size = (0.5 * 72.27/96), 
    linetype = 1, lineend = "butt", arrow = F, inherit.blank = T),
  axis.ticks = element_line(
    colour = "#595A5C", size = (0.5 * 72.27/96), 
    linetype = 1, lineend = "butt", arrow = F, inherit.blank = T)
)

pdf(file.path(out, paste0("samplesOfInterestAll_mergedDensitiesAlt_", format(Sys.Date(), "%Y%m%d"), ".pdf")), 
    useDingbats = FALSE, width = 8, height = 11.5)

par(xpd = NA)

# If you map an aesthetic to a categorical variable, you will get a
# set of contours for each value of that variable

ggplot(data = completeTablesCombined %>% filter(sampleGoodSampleBad %in% c("good", "bad") & !is.na(index))) +
  # geom_density_2d_filled(data = completeTablesCombined %>% dplyr::select(-index, -sampleGoodSampleBad, -TRAIT, -scaledLogLikelihoodRatios),
  #                 aes(y = scaledAltActual, x = scaledPredictions, alpha = after_stat(level)), fill = "grey20", bins = 5) +
  # scale_alpha_ordinal(range = c(0, 0.5)) +
  geom_abline(slope = 1, color="#595A5C", linetype = "solid",
              size = 0.25) +
  geom_abline(intercept = 1, slope = 1, color="grey70", linetype = "dashed",
              size = 0.25) +
  geom_abline(intercept = -1, slope = 1, color="grey70", linetype = "dashed",
              size = 0.25) +
  geom_abline(intercept = 2, slope = 1, color="grey90", linetype = "dashed",
              size = 0.25) +
  geom_abline(intercept = -2, slope = 1, color="grey90", linetype = "dashed",
              size = 0.25) +
  annotate("text",
           label = "+2SD", x = 4.24-2, y = 4, size = 0.5 * 72.27/96 * 6, colour = "grey90", angle = 45) +
  annotate("text",
           label = "+1SD", x = 4.24-1, y = 4, size = 0.5 * 72.27/96 * 6, colour = "grey70", angle = 45) +
  annotate("text",
           label = "-2SD", x = -4.24+2, y = -4, size = 0.5 * 72.27/96 * 6, colour = "grey90", angle = 45) +
  annotate("text",
           label = "-1SD", x = -4.24+1, y = -4, size = 0.5 * 72.27/96 * 6, colour = "grey70", angle = 45) +
  # geom_density_2d(data = completeTablesCombined %>% filter(traitDataType == "continuous") %>% dplyr::select(-index, -sampleGoodSampleBad, -TRAIT),
  #                        aes(y = scaledActual, x = scaledPredictions), alpha = .5, colour = "grey20", size = 0.25) +
  geom_point(data = completeTablesCombined %>% filter(sampleGoodSampleBad %in% c("good", "bad") & !is.na(index)),
             aes(y = scaledAltActual, x = scaledPredictions, colour = sampleGoodSampleBad), shape = 16) +
  geom_segment(aes(y = scaledAltActual, x = scaledPredictions,
                   xend = scaledPredictions, yend = scaledPredictions, colour = sampleGoodSampleBad), alpha = 1, size = 0.25) +
  scale_colour_manual(values = c("good" = "#007E7E", "bad" = "#D13B00"), 
                      name = "Samples", labels = c("good" = "Lowest mix-up predictions", "bad" = "Highest mix-up predictions")) +
  coord_fixed(ratio = 1, xlim = c(-4, 4), ylim = c(-4, 4), expand = TRUE, clip = "on") +
  geom_text_repel(aes(y = scaledAltActual, x = scaledPredictions, label = if_else(traitLabel != "", traitLabel, NA_character_)), na.rm = T,
                  size = 0.5 * 72.27/96 * 8, box.padding = 0.72, max.overlaps = Inf, parse = T,
                  force = 2, force_pull = 1, segment.ncp = 3, segment.angle = 20, colour = "#595A5C") +
  xlab("Predicted phenotype (Scaled and force normalised)") +
  ylab("Predicted phenotype + scaled residuals") +
  ggtitle("Visualizations of scaled residuals reveal foundation of mix-up predictions") +
  facet_grid(~ sampleLabel, labeller = label_parsed) +
  theme(legend.position = "bottom")

dev.off()



pdf(file.path(out, paste0("samplesOfInterestAll_horizontal", format(Sys.Date(), "%Y%m%d"), ".pdf")), 
    useDingbats = FALSE, width = 8, height = 10)

par(xpd = NA)

# ggplot(data = completeTablesCombined,
#        aes(y = scaledActual, x = scaledPredictions, colour = TRAIT)) +
#   geom_density_2d() +
#   facet_grid(rows = vars(traitDataType))

# If you map an aesthetic to a categorical variable, you will get a
# set of contours for each value of that variable
breaks <- c(0, 1, 2, 3, 4, 5)

ggplot(data = completeTablesCombined %>% filter(sampleGoodSampleBad %in% c("bad", "good") & !is.na(index)), aes(x = abs(scaledResiduals), y = TRAIT)) +
  geom_segment(aes(xend = 0, yend = TRAIT), alpha = 1, size = 0.25) +
  geom_point(shape = 16, aes(colour = sampleGoodSampleBad)) +
  # scale_colour_stepsn(colors = RColorBrewer::brewer.pal(5, "Dark2"), breaks = breaks) +
  scale_colour_manual(values = c("good" = "#007E7E", "bad" = "#D13B00"), 
                      name = "Samples", labels = c("good" = "Lowest mix-up predictions", "bad" = "Highest mix-up predictions")) +
  geom_vline(xintercept = 0, color="#595A5C", linetype = "solid",
             size = 0.25) +
  geom_vline(xintercept = 1, color="grey70", linetype = "dashed",
             size = 0.25) +
  geom_vline(xintercept = 2, color="grey90", linetype = "dashed",
             size = 0.25) +
  annotate("text",
           label = "+2SD", vjust = 1, hjust = 0, x = 2.08, y = -0.42, size = 0.5 * 72.27/96 * 6, colour = "grey90") +
  annotate("text",
           label = "+1SD", vjust = 1, hjust = 0, x = 1.08, y = -0.42, size = 0.5 * 72.27/96 * 6, colour = "grey70") +
  coord_cartesian(xlim = c(0, 4), ylim = c(-0.5, NA), expand = TRUE) +
  geom_text_repel(aes(y = TRAIT, x = scaledResiduals, label = if_else(traitLabelSmall != "", traitLabelSmall, NA_character_)), na.rm = T,
                  size = 0.5 * 72.27/96 * 8, parse = T, nudge_y = 0.5 * 72.27/96 * 2, direction = "x",
                  colour = "#595A5C") +
  ylab("Trait") +
  xlab("Scaled residuals") +
  ggtitle("The basis of mix-up predictions") +
  facet_wrap(~sampleGoodSampleBadLabel+indexLabel, labeller = labeller(sampleGoodSampleBadLabel = label_value, indexLabel = label_parsed), ncol = 4) +
  theme(legend.position = "bottom",
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank())

dev.off()


pdf(file.path(out, paste0("samplesOfInterestAll_heatmap", format(Sys.Date(), "%Y%m%d"), ".pdf")), 
    useDingbats = FALSE, width = 8, height = 11.5)

par(xpd = NA)

# If you map an aesthetic to a categorical variable, you will get a
# set of contours for each value of that variable

ggplot(data = completeTablesCombined %>% filter(sampleGoodSampleBad %in% c("good", "bad") & !is.na(index)), aes(x = index, y = TRAIT, fill = abs(scaledResiduals))) +
  geom_tile(colour="white",size=0.2) +
  scale_y_discrete(expand=c(0,0)) +
  scale_x_discrete(expand=c(0,0)) +
  scale_fill_viridis_c(values=c(), option = "C", na.value = "grey90") +
  coord_fixed() +
  ylab("Trait") +
  xlab("Sample") +
  ggtitle("Sample good sample bad?") +
  facet_wrap(~ sampleGoodSampleBad, ncol = 2) +
  theme_classic(base_size=10) +
  theme(legend.position="right",legend.direction="vertical",
        legend.margin=margin(grid::unit(0,"cm")),
        legend.key.height=grid::unit(0.8,"cm"),
        legend.key.width=grid::unit(0.2,"cm"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        plot.margin=margin(0.7,0.4,0.1,0.2,"cm"))
dev.off()
