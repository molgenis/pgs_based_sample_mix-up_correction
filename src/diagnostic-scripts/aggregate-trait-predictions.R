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
library(data.table)
library(stringr)
library(pROC)

##############################
# Define argument parser
##############################

parser <- ArgumentParser(description='')
parser$add_argument('--phenotypes-file',
                    help='path to a tab-delimited file holding all processed phenotype data.')
parser$add_argument('--sample-coupling-file', required = FALSE,
                    help=paste0('file containing genotype sample ids in the first column',
                    'and phenotype sample ids in the second column'))
parser$add_argument('--traits-to-aggregate', required = FALSE,
                    help=paste('file containing trait names to aggregate'))
parser$add_argument('--dir',
                    help='path to directory to recycle results from')
parser$add_argument('--out',
                    help='path to output directory')

##############################
# Define functions
##############################

# Function that returns a matrix of residuals, with phenotype sample across the
# rows, and genotype samples across the columns. 
# This should be able to be read from an existing file:
# - '<intermediateResidualMatrixFileBasePath>/scaledResidualMatrix.tsv', or
# - '<intermediateResidualMatrixFileBasePath>/residualMatrix.rds'.
getLogLikelihoodRatioMatrix <- function(
  intermediateLogLikelihoodRatioMatrixFileBasePath) {
  
  intermediateLogLikelihoodRatioMatrixRdsFilePath <- paste0(
    intermediateLogLikelihoodRatioMatrixFileBasePath, ".logLikelihoodRatios.rds")
  intermediateLogLikelihoodRatioMatrixTsvFilePath <- paste0(
    intermediateLogLikelihoodRatioMatrixFileBasePath, ".logLikelihoodRatios.tsv")
  LogLikelihoodRatioMatrix <- NULL
  
  if (file.exists(intermediateLogLikelihoodRatioMatrixRdsFilePath) 
      && file.access(intermediateLogLikelihoodRatioMatrixRdsFilePath, 4) == 0) {
    
    message(paste0("    Loading log likelihood ratios from '", intermediateLogLikelihoodRatioMatrixRdsFilePath, "'..."))
    LogLikelihoodRatioMatrix <- readRDS(intermediateLogLikelihoodRatioMatrixRdsFilePath)
    
  } else if (file.exists(intermediateLogLikelihoodRatioMatrixTsvFilePath) 
             && file.access(intermediateLogLikelihoodRatioMatrixTsvFilePath, 4) == 0) {
    
    message(paste0("    Loading log likelihood ratios from '", intermediateLogLikelihoodRatioMatrixTsvFilePath, "'..."))
    LogLikelihoodRatioMatrix <- as.matrix(fread(intermediateLogLikelihoodRatioMatrixTsvFilePath), rownames = 1)
    
  } else {
    stop("intermediate log likelihood ratio files not found")
  }
  
  return(LogLikelihoodRatioMatrix)
}

##############################
# Run
##############################
args <- parser$parse_args(c(
  "--phenotypes-file", 
  "/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/data/lifelines/processed/pgs.phenotypes_20201215.ugli.dat",
  "--sample-coupling-file",
  "/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/data/lifelines/processed/pgs.sample-coupling-file.ugli-all.20210218.32817samples.txt",
  "--traits-to-aggregate",
  "/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/jobs/sample-swap-prediction/traits-to-aggregate.txt",
  "--out", 
  "/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/output/sample-swap-prediction/20201120.20210218_ugli-all.NA.80/20210227",
  "--dir", 
  "/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/output/sample-swap-prediction/20201120.20210218_ugli-all.NA.80"))

args <- parser$parse_args(commandArgs(trailingOnly = TRUE))

outputIntermediateStatistics <- F

out <- args$out
dir <- args$dir

traitOutputFilePath <- file.path(dir, "outputStatisticsPerTrait.tsv")

message(strwrap(prefix = " ", initial = "", paste(
  "Loading trait-gwas-mapping:\n", traitOutputFilePath)))

# Load table containing paths for the plink output 
# and corresponding phenotype labels.
traitDescriptionsTable <- fread(
  traitOutputFilePath, 
  quote="", sep = "\t",
  stringsAsFactors=F) %>% 
  select(-V1)

message(strwrap(prefix = " ", initial = "", paste(
  "Loading polygenic scores from:\n", args$trait_gwas_mapping)))

traitsToAggregatePath <- args$traits_to_aggregate
traitDescriptionsTable <- traitDescriptionsTable %>%
  as_tibble() %>%
  inner_join(fread(traitsToAggregatePath, 
                   quote="", header=F, sep = "\t", 
                   col.names=c("trait"),
                   stringsAsFactors=F), by = "trait")

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
                  "Generating sample coupling table from phenotypes file."))
}

predictingInducedMixUps <- F

if (!("original" %in% colnames(link))) {
  link$original <- link$geno
} else {
  predictingInducedMixUps <- T
}

aggregatedLlrMatrix <- matrix(nrow = nrow(link), 
                              ncol = nrow(link), 
                              dimnames = list(link$pheno, link$geno), 0)

aggregatedNumberOfTraits <- matrix(nrow = nrow(link), 
                                   ncol = nrow(link), 
                                   dimnames = list(link$pheno, link$geno), 0)

traitOutputTable <- traitDescriptionsTable

for (traitIndex in 1:nrow(traitDescriptionsTable)) {
  
  polygenicScoreFilePath <- traitDescriptionsTable$polygenicScoreFilePath[traitIndex]
  trait <- traitDescriptionsTable$trait[traitIndex]
  responseDataType <- traitDescriptionsTable$traitDataType[traitIndex]
  traitFileName <- traitDescriptionsTable$traitOutputDir[traitIndex]
  traitDirectory <- file.path(dir, traitFileName)
  
  # Give status update
  message(paste0(traitIndex, " / ", nrow(traitDescriptionsTable), 
                 ": '", trait, "' (", responseDataType, ")."))
  
  naiveBayesParameters <- traitOutputTable[which(traitOutputTable$trait == trait), "naiveBayesParameters"]

  intermediateLogLikelihoodRatioMatrixFileBasePath <- file.path(
    traitDirectory, naiveBayesParameters)
  
  if (!(dir.exists(intermediateLogLikelihoodRatiosFilePath) 
        && traitToInclude)) {
    next
  }
  
  logLikelihoodRatios <- getLogLikelihoodRatioMatrix(
    intermediateLogLikelihoodRatioMatrixFileBasePath)
  
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
  rm(logLikelihoodRatios)
  gc()
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

message(paste0("Exporting output matrix: 'aggregatedLogLikelihoodRatiosDataFrame.tsv'"))

# Export the log likelihood data frame with scaled values
write.table(llrDataFrame, file.path(out, "aggregatedLogLikelihoodRatiosDataFrame.tsv"),
            sep="\t", col.names = T, row.names = F, quote = F)

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
