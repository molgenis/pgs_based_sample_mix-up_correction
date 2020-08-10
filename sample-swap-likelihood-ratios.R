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
library(tidyverse)
library(argparse)
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
parser$add_argument('--sample_coupling_file', required = FALSE,
                    help=paste0('file containing genotype sample ids in the first column',
                    'and phenotype sample ids in the second column'))
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
  
    if (is.factor(covariates[,lab])) {
      newData[, lab] <- sample(covariates[, lab], nrow(newData), replace = T)
    }
    else {
      newData[, lab] <- mean(covariates[,lab])
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
  
  ggsave("sigmoid.pdf", plot = sigmoid)
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

# Function that takes (a set of) estimates, actual values, confounders, 
# and a logistic regression model.

# The function predicts actual values from the independent variables 
# (estimates and confounders) using the given model, 
# and returns the deviance of these predicted values
# to the actual values.
devianceFromLogitRegressionLine <- function(estimate, actual, covariates, logitModel) {
  # Predict actual values using the given model and the supplied independent variables
  predictedActual <- predict(logitModel, 
                             cbind(estimate, covariates), 
                             type = "response")
  
  # Calculate the distance between the actual values and the predicted actual values.
  distanceBetweenSigmoidAndActual <- actual - predictedActual
  
  # Calculate the deviance residuals from the distance between the actual values and
  # the predicted actual values.
  deviance <- sqrt(2 * -log(1 - abs(distanceBetweenSigmoidAndActual)))
  deviance <- sign(distanceBetweenSigmoidAndActual) * deviance
  
  return(unname(deviance))
}


# Function that takes 
residualsFunConstructor <- function(estimate, actual, covariates, responseDataType = "continuous") {
  
  # Return a linear model in case 'actual', is a continuous data type.
  if (responseDataType == "continuous") {
    olsModel <- lm(actual ~ estimate + . + .^2, data = covariates)
    
    print("R-squared:")
    print(summary(olsModel))
    
    # Define residualsFun to use the linear regression model.
    residualsFun <- function(estimate, actual, covariates) {
      return(devianceFromOlsRegressionLine(estimate = estimate, 
                                           actual = actual, 
                                           covariates = covariates, 
                                           olsModel = olsModel))
    }
    
    return(residualsFun)
    
  # Return a logistic model in case 'actual', is a binary data type.
  } else if (responseDataType == "binary") {
    logitModel <- glm(actual ~ estimate + . + .^2, family=binomial(link='logit'), data = covariates)
    
    print(summary(logitModel))
    
    #plotSigmoid(estimate = estimate, actual = actual, covariates = covariates, logitModel = logitModel)
    
    # Define residualsFun to use the logistic regression model.
    residualsFun <- function(estimate, actual, covariates) {
      return(devianceFromLogitRegressionLine(estimate = estimate, 
                                             actual = actual, 
                                             covariates = covariates, 
                                             logitModel = logitModel))
    }
    
    return(residualsFun)
  }
}


# Define function for calculating Z-scores that represents what the Java implementation should do.
calculate.scaledResiduals <- function(estimate, actual, covariates, sampleNames = NULL, responseDataType = "continuous") {
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
  
  # Create names
  if (is.null(sampleNames)) {
    sampleNames <- 1:length(actual)
  }
  
  actualDataFrame <- as.data.frame(c(list(phenotypeSamples = sampleNames, actual = actual), covariates))
  estimateDataFrame <- data.frame(genotypeSamples = sampleNames, estimate = estimate)
  scaledResidualDataFrame <- expand_grid(actualDataFrame, estimateDataFrame)
  scaledResidualDataFrame$scaledResiduals <- residualsFun(estimate = scaledResidualDataFrame$estimate, 
                                                         actual = scaledResidualDataFrame$actual, 
                                                         covariates = scaledResidualDataFrame[colnames(covariates)])
  
  scaledResidualDataFrame$scaledResiduals <- (scaledResidualDataFrame$scaledResiduals - residuals.mean) / residuals.sd
  
  return(
    scaledResidualDataFrame %>% 
    select(phenotypeSamples, genotypeSamples, scaledResiduals))
}

# Function for converting a matrix of scaled residuals to log likelihood ratios
scaledResidualsToLlr.naiveBayes <- function(scaledResiduals, group, nBins = 50) {

  # Extract the null-residuals; 
  # the residuals belonging to the matches that are assumed to be correct.
  nullResiduals <- scaledResiduals[group == "null"]
  
  # Extract the alternative residuals; 
  # the residuals belonging to the matches that are assumed to be sample-swaps.
  alternativeResiduals <- scaledResiduals[group == "alternative"]
  
  # nullDensity <- density(
  #   nullResiduals, 
  #   bw = "SJ",
  #   from = 0, to = max(scaledResiduals))
  # alternativeDensity <- density(
  #   alternativeResiduals, 
  #   bw = "SJ",
  #   from = 0, to = max(scaledResiduals))
  # 
  # plot(nullDensity)
  # lines(alternativeDensity)
  
  nullTiles <- ntile(nullResiduals, nBins)
  breaks <- sapply(1:nBins, function(bin) min(nullResiduals[nullTiles == bin]))
  
  breaks[1] <- min(scaledResiduals - 1)
  breaks <- c(breaks, max(scaledResiduals) + 1)
  
  alternativeTiles <- cut(alternativeResiduals, breaks = breaks, labels = FALSE)
  
  nullLikelihoods <- sapply(
    1:nBins, 
    function(bin) sum(nullTiles == bin) / length(nullTiles))
  
  alternativeLikelihoods <- sapply(
    1:nBins, 
    function(bin) sum(alternativeTiles == bin) / length(alternativeTiles))
  
  likelihoodRatioMap <- alternativeLikelihoods / nullLikelihoods
  
  allTiles <- cut(scaledResiduals, breaks = breaks, labels = FALSE)
  likelihoodRatios <- sapply(allTiles, function(tile) likelihoodRatioMap[tile])
  
  return(log(likelihoodRatios))
}

# Function for converting a matrix of scaled residuals to likelihood ratios.
# This function employs a gaussian naive bayes method to calculate likelihoods.
scaledResidualsToLlr.gaussianNaiveBayes <- function(scaledResiduals, group, error = 0.01) {
  
  # Extract the null-residuals; 
  # the residuals belonging to the matches that are assumed to be correct.
  nullResiduals <- scaledResiduals[group == "null"]
  
  # Extract the alternative residuals; 
  # the residuals belonging to the matches that are assumed to be sample-swaps.
  alternativeResiduals <- scaledResiduals[group == "alternative"]
  
  # Calculate the mean of both groups of residuals.
  nullMean <- mean(nullResiduals)
  alternativeMean <- mean(alternativeResiduals)
  
  # Calculate the standard deviation in both groups of residuals.
  nullSd <- sd(nullResiduals)
  alternativeSd <- sd(alternativeResiduals)
  
  # Declare function for calculating the difference between two log-transformed
  # likelihoods. returns -Inf if l1 and l2 are (approximately) equal. l1 must be > l2.
  # retrieved from: https://stats.stackexchange.com/a/383524 on 22-07-2020.
  subtractLogTransformedProbabilities <- function(l1, l2) {
    return(l1 + log1p(-exp(-(l1 - l2))))
  }
  
  # Declare function for calculating the log likelihood of 'value' being sampled from 
  # a normal distribution with the given mean and standard deviation 'sd'.
  gaussianLogLikelihood <- function(value, mean, sd) {
    return(subtractLogTransformedProbabilities(
      pnorm(value + (error / 2), mean = mean, sd = sd, log.p = TRUE),
      pnorm(value - (error / 2), mean = mean, sd = sd, log.p = TRUE)
      ))
  }
  
  # Calculate, for every residual, the likelihood of the residual being sampled from
  # the corresponding normal distribution of null residuals.
  nullLogLikelihoods <- sapply(
    nullResiduals, gaussianLogLikelihood, nullMean, nullSd)
  
  # Calculate, for every residual, the likelihood of the residual being sampled from
  # the corresponding normal distribution of alternative residuals.
  alternativeLogLikelihoods <- sapply(
    alternativeResiduals, gaussianLogLikelihood, alternativeMean, alternativeSd)
  
  # For every residual, calculate the likelihood ratio the residual belonging to a
  # sample swap.
  logLikelihoodRatios <- alternativeLogLikelihoods - nullLogLikelihoods
  
  return(logLikelihoodRatios)
}

# Define function for calculating the AUC
calculate.auc <- function(actual, predictor) {
  pr <- prediction(predictor, actual)
  
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

##############################
# Run
##############################
args <- parser$parse_args(commandArgs(trailingOnly = TRUE))
# args <- parser$parse_args(c("--profiles",
#                          "~/pgs_based_mixup_correction/jobs/pgs_output_processing/mapping_with_responseType.txt",
#                          "--phenotypes_path",
#                          "/home/umcg-rwarmerdam/pgs_based_mixup_correction/data/lldeep/samples/LLDeep_GoNL_samples_V01_20200313.txt"))

# Load table containing paths for the plink output 
# and corresponding phenotype labels.
traitDescriptionsTable <- read.table(
  args$trait_gwas_mapping, quote="", header=T, sep = "\t",
  col.names=c("trait", "traitDataType", "summaryStatistics"), stringsAsFactors=F)

# Get the paths to the polygenic scores.
basePathWithPolygenicScores <- args$base_pgs_path
traitDescriptionsTable$polygenicScoreFilePath <- file.path(basePathWithPolygenicScores, 
                                                 traitDescriptionsTable$summaryStatistics, 
                                                 "full.UGLI.pgs.profile")

# Get the plink output paths
polygenicScoreFilePaths <- traitDescriptionsTable$polygenicScoreFilePath

# Get the phenotype labels
traits <- traitDescriptionsTable$trait

# Get the output path
out <- args$out

# Load the phenotypes 
phenotypesFilePath <- args$phenotypes_file
phenotypesTable <- read.table(phenotypesFilePath, header=T, quote="", sep="\t",
                              col.names = c("ID", "AGE", "SEX", "VALUE", "TRAIT"))

link <- data.frame(geno = phenotypesTable$ID, pheno = phenotypesTable$ID)

# Get the link path
if (!is.null(args$sample_coupling_file)) {
  link <- fread(args$sample_coupling_file, stringsAsFactors=F,
                col.names = c("geno", "pheno"))
}

llRdataFrameList <- list()
pearson.correlations <- data.frame(trait = traits,
                                   pearson.not_corrected = 0.0, 
                                   pearson.corrected.sex = 0.0,
                                   pearson.corrected.age = 0.0,
                                   pearson.corrected.both = 0.0)
rownames(pearson.correlations) <- pearson.correlations$phenotype

# Loop trough traits

for (fileIndex in c(1:length(polygenicScoreFilePaths))) {
  
  polygenicScoreFilePath <- polygenicScoreFilePaths[fileIndex]
  trait <- traits[fileIndex]
  responseDataType <- traitDescriptionsTable$traitDataType[fileIndex]
  
  phenotypeTable <- phenotypesTable %>%
    filter(TRAIT == trait) %>%
    rename(pheno = ID) %>%
    inner_join(link, by="pheno")
  
  # Give status update
  message(paste0(fileIndex, " / ", length(profilePaths), ": '", trait, "' (", responseDataType, ")."))
  
  if (responseDataType == "binary") {
    phenotypeFrequencyTable <- table(phenotypeTable)
    
    message(paste0("    Available for ", nrow(phenotypeTable), 
                   " samples (number of 0's = ", phenotypeFrequencyTable["0"],
                   ", 1's = ", phenotypeFrequencyTable["1"], ")."))
    warning("Skipping...")
    next
  } else if (responseDataType == "continuous") {
    message(paste0("    Available for ", nrow(phenotypeTable), "samples."))
  }

  message(paste0("    Loading polygenic scores from '", polygenicScoreFilePaths, "'..."))
  
  # Read the PLINK polygenic score table.
  polygenicScores <- read.table(
    polygenicScoreFilePaths,
    header=T) %>%
    rename(SCORESUM = "PGS")
  
  completeTable <- phenotypeTable %>%
    inner_join(polygenicScores, by = c("geno" = "IID"))
  
  print(head(completeTable))
  
  initial.cor.test.results <- cor.test(completeTable$VALUE, completeTable$PGS)
  print(paste0("Initial R-squared of correlation = ", initial.cor.test.results$estimate ^ 2))
  pearson.correlations[fileIndex, "pearson.not_corrected"] <- initial.cor.test.results$estimate
  
  # Correct for age and sex
  trait2pgs.corrected <- completeTable
  model.complete <- lm(VALUE ~ AGE + SEX + AGE * SEX, data = trait2pgs.corrected)
  trait2pgs.corrected$actual <- resid(model.complete)
  
  # Scale both the actual and estimated traits
  #trait2pgs.corrected$actual <- rank(trait2pgs.corrected$actual)/length(trait2pgs.corrected$actual)
  #trait2pgs.corrected$PGS <- rank(trait2pgs.corrected$PGS)/length(trait2pgs.corrected$PGS)
  #trait2pgs.corrected$actual <- scale(trait2pgs.corrected$actual)
  #trait2pgs.corrected$PGS <- scale(trait2pgs.corrected$PGS)
  
  # Output the correlation of the corrected traits
  corrected.cor.test.results <- cor.test(trait2pgs.corrected$actual, trait2pgs.corrected$PGS)
  print(paste0("R-squared of corrected traits = ", corrected.cor.test.results$estimate ^ 2))
  pearson.correlations[fileIndex, "pearson.corrected.both"] <- corrected.cor.test.results$estimate
  
  # Calculate z-score matrix based on polygenic scores and actual phenotypes,
  # using the chosen function for calculating residuals.
  scaledResidualsDataFrame <- calculate.scaledResiduals(
    estimate = completeTable$PGS, 
    actual = completeTable$VALUE, 
    covariates = completeTable[c("AGE", "SEX")],
    sampleNames = completeTable$geno,
    responseDataType = responseDataType)
  
  scaledResidualsDataFrame$group <- "alternative"
  
  scaledResidualsDataFrame$group[which(
    scaledResidualsDataFrame$genotypeSamples == scaledResidualsDataFrame$phenotypeSamples)] <- "null"
  
  scaledResidualsDataFrame$group <- factor(scaledResidualsDataFrame$group, c("null", "alternative"))
  
  # Calculate the likelihood ratios for every residual being from the distribution of possible mix-ups.
  
  scaledResidualsDataFrame$logLikelihoodRatios <- scaledResidualsToLlr.naiveBayes(
    scaledResidualsDataFrame$scaledResiduals,
    group = scaledResidualsDataFrame$group,
    nBins = 128)
  
  # scaledResidualsDataFrame$logLikelihoodRatios <- scaledResidualsToLlr.gaussianNaiveBayes(
  #   scaledResidualsDataFrame$scaledResiduals,
  #   group = scaledResidualsDataFrame$group)
  
  llRdataFrameList[[fileIndex]] <- scaledResidualsDataFrame %>% 
    select(genotypeSamples, phenotypeSamples, group, logLikelihoodRatios)
  
  ggplot(scaledResidualsDataFrame, aes(x=scaledResiduals, stat(density), fill=group)) +
    geom_histogram(bins = 72, alpha=.5, position="identity") +
    geom_line(aes(y=logLikelihoodRatios)) +
    xlab("Scaled residuals") + ggtitle(paste0("Scaled residuals for trait '", trait, "'"))
  
  ggsave(file.path(out, trait, "/scaledResidualsHistogram.png"), width=8, height=7)
}

mergedLlrDataFrame <- 
  Reduce(function(...) merge(..., by = c("genotypeSamples", "phenotypeSamples", "group")), llRdataFrameList)

lRProducts <- data.frame(mergedLlrDataFrame[1:3], logLikelihoodRatios = apply(mergedLlrDataFrame[-1:-3], 1, sum))

message(paste0("Calculated overall AUC: ", calculate.auc(lRProducts$group == "null", lRProducts$logLikelihoodRatios)))

# zscore.sums$zscore.avg <- zscore.sums$zscore.sum / sqrt(length(zscore.list))
# zscore.sums$zscore.avg <- zscore.sums$zscore.sum

# Write these Z-scores for later.
write.table(lRProducts, file.path(out, "/likelihoodRatios.tsv"),
            sep="\t", col.names = T, row.names = T, quote = F)

ggplot(lRProducts, aes(x=logLikelihoodRatios, stat(density), fill=group)) +
  geom_histogram(bins = 32, alpha=.5, position="identity") +
  xlab("Likelihood ratios") + ggtitle(paste0("LR overall"))

ggsave(file.path(out, "/likelihoodRatioHistogram.png"), width=8, height=7)

ggplot(lRProducts, aes(x=group, y=logLikelihoodRatios)) +
  geom_boxplot() + ggtitle(paste0("Likelihood ratio distributions comparison overall"))

ggsave(file.path(out, "/likelihoodRatioBoxplot.png"), width=8, height=7)

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
