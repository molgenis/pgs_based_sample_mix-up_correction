#!/usr/bin/env Rscript

#############################################################
## c.a.warmerdam@umcg.nl - April 2020
## Processes plink polygenic scores and phenotypes and
## compares these between samples to quantify the
## to what degree samples are swapped with each other.
#############################################################

##############################
# Load librairies
##############################
library(tidyverse)
library(argparse)
library(ROCR)

##############################
# Define argument parser
##############################

parser <- ArgumentParser(description='')
parser$add_argument('--profiles',
                    help='path to plink --score output')
parser$add_argument('--phenotypes_path',
                    help='phenotype table')
parser$add_argument('--sample_coupling_file', required = FALSE,
                    help=paste0('file containing genotype sample ids in the first column',
                    'and phenotype sample ids in the second column'))

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

# Function for converting a matrix of scaled residuals to likelihood ratios
scaledResidualsToLr <- function(scaledResiduals, group, nBins = 50) {

  nullResiduals <- scaledResiduals[group == "null"]
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
  
  return(likelihoodRatios)
}

# Define function for calculating the AUC
calculate.auc <- function(actual, predictor) {
  model <- glm(actual ~ predictor,
               family=binomial(link='logit'))
  
  p <- predict(model, type="response")
  pr <- prediction(p, actual)
  
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
traitDescriptionsTable <- read.csv(
  args$profiles, quote="", header=F, 
  col.names=c("file", "phenotype", "traitDataType"), stringsAsFactors=F)

# Get the plink output paths
profiles <- traitDescriptionsTable$file

# Get the phenotype labels
phenotypes <- traitDescriptionsTable$phenotype

# Load the phenotypes 
phenotypes_path <- args$phenotypes_path
phenotypes_table <- read.table(phenotypes_path, header=T, row.names=1, quote="", sep="\t")

link <- data.frame(geno = rownames(phenotypes_table), pheno = rownames(phenotypes_table))

# Get the link path
if (!is.null(args$sample_coupling_file)) {
  link <- fread(args$sample_coupling_file, stringsAsFactors=F,
                col.names = c("geno", "pheno"))
}

# Get the tables for correcting phenotypes
trait.table <- phenotypes_table[,"geslacht", drop=F]
trait.name <- colnames(trait.table)[1]
colnames(trait.table) <- c("crrct.trait")

age.table <- phenotypes_table[,"age_bl1", drop=F]

lRdataFrameList <- list()
pearson.correlations <- data.frame(phenotype = phenotypes,
                                   pearson.not_corrected = 0.0, 
                                   pearson.corrected.sex = 0.0,
                                   pearson.corrected.age = 0.0,
                                   pearson.corrected.both = 0.0)
rownames(pearson.correlations) <- pearson.correlations$phenotype

# Loop trough traits

for (file.index in c(1:length(profiles))) {
  
  filename <- profiles[file.index]
  print(filename)
  print(phenotypes[file.index])
  all.trait <- phenotypes_table[,phenotypes[file.index], drop=F]
  colnames(all.trait) <- c("actual")
  
  all.trait <- all.trait %>% 
    mutate(pheno = rownames(all.trait)) %>%
    inner_join(link, by="pheno") %>%
    select(geno, actual)
  
  rownames(all.trait) <- all.trait$geno
  
  responseDataType <- traitDescriptionsTable$traitDataType[file.index]

  # Read the PLINK polygenic score table.
  all.PGSs <- read.table(
    filename,
    header=T, row.names="IID")
  
  mat.all.PGSs <- merge(all.trait, all.PGSs, by="row.names")
  mat.all.PGSs.colnames <- colnames(mat.all.PGSs)
  mat.all.PGSs.colnames[1] <- "ID"
  mat.all.PGSs.colnames[mat.all.PGSs.colnames=="SCORESUM"] <- "PGS"
  
  colnames(mat.all.PGSs) <- mat.all.PGSs.colnames
  rownames(mat.all.PGSs) <- mat.all.PGSs$ID
  
  head(mat.all.PGSs)
  
  trait2pgs <- mat.all.PGSs[c("actual", "PGS")]
  
  print(paste0("Processing '", tools::file_path_sans_ext(basename(filename)), "'..."))
  colnames(trait2pgs) <- c("actual", "PGS")
  
  initial.cor.test.results <- cor.test(trait2pgs$actual, trait2pgs$PGS)
  print(paste0("Initial R-squared of correlation = ", initial.cor.test.results$estimate ^ 2))
  pearson.correlations[file.index, "pearson.not_corrected"] <- initial.cor.test.results$estimate
  
  # Merge the table of polygenic scores vs actual traits, and the table with a trait to correct for
  combined <- merge(trait2pgs, trait.table, by="row.names")
  rownames(combined) <- combined$Row.names
  combined$Row.names <- NULL
  combined <- merge(combined, age.table, by="row.names")
  
  combined <- combined[!is.na(combined$actual),]
  
  # Perform actual correction
  # trait2pgs.corrected <- combined %>% group_by(crrct.trait) %>% mutate(actual = scale(actual))
  
  # Correct for sex
  # trait2pgs.sexCorrected <- combined
  # model.geslacht <- lm(actual ~ crrct.trait, data = combined)
  # trait2pgs.sexCorrected$actual <- resid(model.geslacht)
  # 
  # sex.cor.test.results <- cor.test(trait2pgs.sexCorrected$actual, trait2pgs.sexCorrected$PGS)
  # print(paste0("Sex corrected pearson correlation = ", sex.cor.test.results$estimate))
  # pearson.correlations[file.index, "pearson.corrected.sex"] <- sex.cor.test.results$estimate
  
  # Correct for age
  # trait2pgs.ageCorrected <- combined
  # model.age <- lm(actual ~ age_bl1, data = combined)
  # trait2pgs.ageCorrected$actual <- resid(model.age)
  
  # age.cor.test.results <- cor.test(combined$actual, combined$PGS)
  # print(paste0("Age corrected pearson correlation = ", age.cor.test.results$estimate))
  # pearson.correlations[file.index, "pearson.corrected.age"] <- age.cor.test.results$estimate
  
  # Correct for age and sex
  trait2pgs.corrected <- combined
  model.complete <- lm(actual ~ age_bl1 + crrct.trait + age_bl1 * crrct.trait, data = combined)
  trait2pgs.corrected$actual <- resid(model.complete)
  
  # Scale both the actual and estimated traits
  #trait2pgs.corrected$actual <- rank(trait2pgs.corrected$actual)/length(trait2pgs.corrected$actual)
  #trait2pgs.corrected$PGS <- rank(trait2pgs.corrected$PGS)/length(trait2pgs.corrected$PGS)
  #trait2pgs.corrected$actual <- scale(trait2pgs.corrected$actual)
  #trait2pgs.corrected$PGS <- scale(trait2pgs.corrected$PGS)
  
  # Output the correlation of the corrected traits
  corrected.cor.test.results <- cor.test(trait2pgs.corrected$actual, trait2pgs.corrected$PGS)
  print(paste0("R-squared of corrected traits = ", corrected.cor.test.results$estimate ^ 2))
  pearson.correlations[file.index, "pearson.corrected.both"] <- corrected.cor.test.results$estimate
  
  # Calculate zscore matrix based on polygenic scores and actual phenotypes,
  # using the chosen function for calculating residuals.
  scaledResidualsDataFrame <- calculate.scaledResiduals(
    estimate = combined$PGS, 
    actual = combined$actual, 
    covariates = combined[c("age_bl1", "crrct.trait")],
    sampleNames = combined$Row.names,
    responseDataType = responseDataType)
  
  scaledResidualsDataFrame$group <- "alternative"
  
  scaledResidualsDataFrame$group[which(
    scaledResidualsDataFrame$genotypeSamples == scaledResidualsDataFrame$phenotypeSamples)] <- "null"
  
  scaledResidualsDataFrame$group <- factor(scaledResidualsDataFrame$group, c("null", "alternative"))
  
  scaledResidualsDataFrame$likelihoodRatios <- scaledResidualsToLr(
    scaledResidualsDataFrame$scaledResiduals,
    group = scaledResidualsDataFrame$group,
    nBins = 128)
  
  lRdataFrameList[[file.index]] <- scaledResidualsDataFrame %>% 
    select(genotypeSamples, phenotypeSamples, group, likelihoodRatios)
  
  ggplot(scaledResidualsDataFrame, aes(x=scaledResiduals, stat(density), fill=group)) +
    geom_histogram(bins = 72, alpha=.5, position="identity") +
    geom_line(aes(y=likelihoodRatios)) +
    xlab("Scaled residuals") + ggtitle(paste0("Scaled residuals for trait '", phenotypes[file.index], "'"))
  
  ggsave(paste0(dirname(filename), "/scaledResidualsHistogram.png"), width=8, height=7)
}

basedir <- dirname(dirname(filename))

# merged.zscore.df <- 
#   Reduce(function(...) merge(..., by = c("Var1", "Var2", "testedMatch")), zscore.list)

# print(str(merged.zscore.df))
# zscore.sums <- data.frame(merged.zscore.df[1:3], zscore.sum = apply(merged.zscore.df[-1:-3], 1, sum))

mergedLRdataFrame <- 
  Reduce(function(...) merge(..., by = c("genotypeSamples", "phenotypeSamples", "group")), lRdataFrameList)

lRProducts <- data.frame(mergedLRdataFrame[1:3], likelihoodRatios = apply(mergedLRdataFrame[-1:-3], 1, prod))

print(paste0("Calculated overall AUC: ", calculate.auc(lRProducts$group == "null", lRProducts$likelihoodRatios)))

# zscore.sums$zscore.avg <- zscore.sums$zscore.sum / sqrt(length(zscore.list))
# zscore.sums$zscore.avg <- zscore.sums$zscore.sum

# Write these Z-scores for later.
write.table(lRProducts, paste0(basedir, "/likelihoodRatios.tsv"),
            sep="\t", col.names = T, row.names = T, quote = F)

ggplot(lRProducts, aes(x=likelihoodRatios, stat(density), fill=group)) +
  geom_histogram(bins = 32, alpha=.5, position="identity") +
  xlab("Likelihood ratios") + ggtitle(paste0("LR overall"))

ggsave(paste0(basedir, "/likelihoodRatioHistogram.png"), width=8, height=7)

ggplot(lRProducts, aes(x=group, y=likelihoodRatios)) +
  geom_boxplot() + ggtitle(paste0("Likelihood ratio distributions comparison overall"))

ggsave(paste0(basedir, "/likelihoodRatioBoxplot.png"), width=8, height=7)

# Pearson correlations
# pearson.correlations <- pearson.correlations[order(pearson.correlations$pearson.corrected),]
pearson.correlations$r2.not_corrected <- pearson.correlations$pearson.not_corrected ^ 2
pearson.correlations$r2.corrected.sex <- pearson.correlations$pearson.corrected.sex ^ 2
pearson.correlations$r2.corrected.age <- pearson.correlations$pearson.corrected.age ^ 2
pearson.correlations$r2.corrected.both <- pearson.correlations$pearson.corrected.both ^ 2

write.table(pearson.correlations, paste0(dirname(filename), "/correlations_comparison.tsv"))

pearson.correlations.melted <- melt(pearson.correlations[,c('phenotype','r2.not_corrected','r2.corrected.sex', 'r2.corrected.age', 'r2.corrected.both')], value.name = "r2")
pearson.correlations.melted <- pearson.correlations.melted[order(pearson.correlations.melted[,"r2"]),]

ggplot(pearson.correlations.melted,aes(x = phenotype,y = r2)) +
  geom_bar(aes(fill = variable),stat = "identity",position = "dodge") +
  coord_flip()

ggsave(paste0(dirname(filename), "/correlations_comparison.png"), width=8, height=7)

# Calculate the correlation of residual Z-scores between phenotypes
res.zscore.correlations <- sapply(lRdataFrameList, function(x) sapply(lRdataFrameList, function(y) cor(x["likelihoodRatios"], y["likelihoodRatios"])))
rownames(res.zscore.correlations) <- phenotypes
colnames(res.zscore.correlations) <- phenotypes
res.zscore.correlations <- as.matrix(res.zscore.correlations)

res.zscore.correlations.melted <- melt(res.zscore.correlations, value.name="pearson.cor")

# Color Brewer palette
ggplot(data = res.zscore.correlations.melted, aes(x=Var1, Var2, fill=pearson.cor)) +
  ggtitle("Correlations of LRs between phenotypes") +
  xlab("Phenotypes") +
  ylab("Phenotypes") +
  geom_tile() +
  scale_fill_distiller(palette = "Blues") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(paste0(dirname(filename), "/residual_LikelihoodRatioCorrelations.png"), width=8, height=7)
