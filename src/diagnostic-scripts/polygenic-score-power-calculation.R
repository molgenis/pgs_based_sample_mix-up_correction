#!/usr/bin/env Rscript

#############################################################
## c.a.warmerdam@umcg.nl - November 2020
## Processes plink polygenic scores and phenotypes and
## calculates predictive power for every trait
#############################################################

##############################
# Load libraries
##############################
library(tidyverse)
library(argparse)
library(data.table)
library(pROC)

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

# Define function for calculating the AUC
calculate.auc <- function(actual, predictor) {
  pr <- prediction(predictor, actual)
  
  performance <- performance(pr, measure = "auc")
  auc <- performance@y.values[[1]]
  return(auc)
}

##############################
# Run
##############################
# args <- parser$parse_args(c("--trait-gwas-mapping", "/groups/umcg-lld/tmp01/other-users/umcg-rwarmerdam/pgs_based_mixup_correction/scripts/r-scripts/pgs_based_sample_mix-up_correction/trait-gwas-mapping.txt",
#                             "--sample-coupling-file", "/home/umcg-rwarmerdam/pgs_based_mixup_correction-ugli/data/lifelines/processed/pgs.sample-coupling-file.ugli.20201014.perm_5120samples_51mixUps.txt",
#                             "--base-pgs-path", "/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/output/PRScs/20201120/",
#                             "--phenotypes-file", "/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/data/lifelines/processed/pgs.phenotypes.ugli.dat",
#                             "--out", "/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/output/sample-swap-prediction/20200811.test/"))
args <- parser$parse_args(c("--trait-gwas-mapping", "/groups/umcg-lld/tmp01/other-users/umcg-rwarmerdam/pgs_based_mixup_correction/scripts/r-scripts/pgs_based_sample_mix-up_correction/trait-gwas-mapping.txt",
                            "--base-pgs-path", "/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/output/PRScs/20201120/",
                            "--phenotypes-file", "/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/data/lifelines/processed/pgs.phenotypes_20201215.ugli.dat",
                            "--out", "/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/output/polygenic-score-power-calculation/20201120.20201218/"))

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

if (!dir.create(out, recursive = T)) {
  warning(paste0("Could not create directory '", out, "'"))
}

message(strwrap(prefix = " ", initial = "", paste(
  "Loading phenotype table:\n", args$phenotypes_file)))

# Load the phenotypes 
phenotypesFilePath <- args$phenotypes_file
phenotypesTable <- fread(phenotypesFilePath, header=T, quote="", sep="\t") %>%
  distinct() %>%
  rename_all(recode, "UGLI_ID" = "ID") %>%
  mutate(ID = as.character(ID),
         SEX = case_when(SEX == 1 ~ "Female", SEX == 2 ~ "Male", TRUE ~ as.character(SEX)),
         SEX = factor(SEX, levels = c("Female", "Male"))) %>%
  group_by(ID) %>%
  filter(!(any(AGE < 18)) & !is.na(VALUE)) %>%
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

traitDescriptionsTable <- traitDescriptionsTable %>%
  mutate(pearsonCorrelation = NA_real_,
         spearmanCorrelation = NA_real_,
         rocCurveAuc = NA_real_,
         pearsonCorrelationOnAdjustedPhenotype = NA_real_,
         spearmanCorrelationOnAdjustedPhenotype = NA_real_,
         rocCurveAucOnAdjustedPhenotype = NA_real_)

# Loop trough traits
for (traitIndex in 1:nrow(traitDescriptionsTable)) {
  
  polygenicScoreFilePath <- traitDescriptionsTable$polygenicScoreFilePath[traitIndex]
  trait <- traitDescriptionsTable$trait[traitIndex]
  responseDataType <- traitDescriptionsTable$traitDataType[traitIndex]
  
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
  
  traitDescriptionsTable[traitIndex, "pearsonCorrelation"] <- 
    cor.test(completeTable$PGS, completeTable$VALUE, method = "pearson")$estimate
  traitDescriptionsTable[traitIndex, "spearmanCorrelation"] <- 
    cor.test(completeTable$PGS, completeTable$VALUE, method = "spearman")$estimate
  if (length(unique(completeTable$VALUE)) == 2) {
    traitDescriptionsTable[traitIndex, "rocCurveAuc"] <- 
      auc(completeTable$VALUE, completeTable$PGS)
  }
  
  correctionModel <- lm(VALUE ~ AGE + SEX + AGE:SEX, data = completeTable)
  completeTable$VALUE.CORRECTED <- resid(correctionModel)
  
  traitDescriptionsTable[traitIndex, "pearsonCorrelationOnAdjustedPhenotype"] <- 
    cor.test(completeTable$PGS, completeTable$VALUE.CORRECTED, method = "pearson")$estimate
  traitDescriptionsTable[traitIndex, "spearmanCorrelationOnAdjustedPhenotype"] <- 
    cor.test(completeTable$PGS, completeTable$VALUE.CORRECTED, method = "spearman")$estimate
  if (length(unique(completeTable$VALUE.CORRECTED)) == 2) {
    traitDescriptionsTable[traitIndex, "rocCurveAucOnAdjustedPhenotype"] <- 
      auc(completeTable$VALUE, completeTable$PGS)
  }
}

write.table(traitDescriptionsTable, file.path(out, "outputStatisticsPerTrait.tsv"),
            sep="\t", col.names = T, row.names = T, quote = F)