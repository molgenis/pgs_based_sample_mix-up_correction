#!/usr/bin/env Rscript

#############################################################
## c.a.warmerdam@umcg.nl - January 2021
## Apply threshold and filter on sex correspondence
#############################################################

##############################
# Load libraries
##############################
library(tidyverse)
library(argparse)
library(data.table)
library(pROC)

##############################
# Define ggplot theme
##############################

old <- theme_set(theme_classic())
theme_update(line = element_line(
  colour = "black", size = (0.5 * ggplot2::.pt * 72.27/96), 
  linetype = 1, lineend = "butt", arrow = F, inherit.blank = T),
  strip.background = element_rect(colour = NA, fill = NA),
  axis.line = element_line(
    colour = "#595A5C", size = (0.5 * ggplot2::.pt * 72.27/96), 
    linetype = 1, lineend = "butt", arrow = F, inherit.blank = T),
  axis.ticks = element_line(
    colour = "#595A5C", size = (0.5 * ggplot2::.pt * 72.27/96), 
    linetype = 1, lineend = "butt", arrow = F, inherit.blank = T)
)

##############################
# Define argument parser
##############################

parser <- ArgumentParser(description='')
parser$add_argument('--dir',
                    help='path to directory to recycle results from')

##############################
# Run
##############################

# Parse arguments.
args <- parser$parse_args(commandArgs(trailingOnly = TRUE))

dir <- args$dir

idefixPredictionPath <- file.path(dir, "idefixPredictions.txt")

if (file.exists(idefixPredictionPath)) {
  
  # Load the idefix predictions
  results <- read_delim(idefixPredictionPath, delim = "\t")
  
} else {
  
  logLikelihoodMatrix <- as.matrix(fread(file.path(dir, "aggregatedLogLikelihoodRatiosMatrix.tsv")), rownames = 1)
  numberOfTraitMatrix <- as.matrix(fread(file.path(dir, "aggregatedNumberOfTraitsMatrix.tsv")), rownames = 1)
  
  scaledLogLikelihoodMatrix <- t(apply(logLikelihoodMatrix, 1, function(x) scale(x)))
  
  # Extract the diagonal values
  diagValues <- scaledLogLikelihoodMatrix[lower.tri(scaledLogLikelihoodMatrix, diag = TRUE)
                                          & upper.tri(scaledLogLikelihoodMatrix, diag = TRUE)]
  
  # Extract the diagonal values
  diagTraitNumbers <- numberOfTraitMatrix[lower.tri(numberOfTraitMatrix, diag = TRUE)
                                          & upper.tri(numberOfTraitMatrix, diag = TRUE)]
  
  # Extract the alternative residuals; 
  # the residuals belonging to the matches that are assumed to be sample-swaps.
  offDiagValues <- scaledLogLikelihoodMatrix[lower.tri(scaledLogLikelihoodMatrix) 
                                             | upper.tri(scaledLogLikelihoodMatrix)]
  
  results <- tibble(scaledLogLikelihoodRatios = diagValues, 
                    numberOfTraits = diagTraitNumbers, 
                    ID = rownames(scaledLogLikelihoodMatrix))
  
  
  write.table(results, "idefix_predictions.txt", row.names = F, col.names = T, quote = F, sep = "\t")
}

# 
link <- read.table("/groups/umcg-lifelines/tmp01/releases/gsa_linkage_files/v1/gsa_linkage_file.dat", header = T)
results[results$scaledLogLikelihoodRatios > 3.25157236407722,]
genotyping <- genotyping %>% inner_join(link, by = "UGLI_ID")

qc <- fread("/groups/umcg-lifelines/tmp01/releases/gsa_genotypes/v1/Logs/ugli_qc_release1_samples.csv")

plate_layout <- str_match(qc$Sample_ID, regex_string)
colnames(plate_layout) <- c("Sample_ID", "Plate", "Ltr", "Nmbr")
plate_layout <- as_tibble(plate_layout)

qc <- qc %>% inner_join(plate_layout, by = "Sample_ID")


qc_full <- genotyping %>% inner_join(qc, by = c("genotyping_name" = "Sample_ID"))

resultsWithPcaEur <- results %>%
  inner_join(gsaLink, by = "UGLI_ID") %>%
  inner_join(qc %>% select(Sample_ID, PCA_european), by = c("genotyping_name" = "Sample_ID")) %>%
  mutate(isOfEuropeanDescent = factor(case_when(PCA_european ~ "European", 
                                                !PCA_european ~ "non-European"),
                                      levels = c("European", "non-European")))

resultsEurOnly <- resultsWithPcaEur %>% filter(PCA_european)
nrow(resultsEurOnly %>% filter(scaledLogLikelihoodRatios > 2.612))

t.test(scaledLogLikelihoodRatios ~ isOfEuropeanDescent, data = resultsWithPcaEur)

pdf(file.path(out, paste0("thresholding_", format(Sys.Date(), "%Y%m%d"), ".pdf")), 
    useDingbats = FALSE, width = 8, height = 4)

par(xpd = NA)

ggplot(resultsWithPcaEur, aes(x = scaledLogLikelihoodRatios, fill = isOfEuropeanDescent)) +
  geom_histogram(binwidth = 0.4) +
  scale_fill_manual(values = c("#119696", "#F95B1C")) +
  geom_rug(data = resultsWithPcaEur) +
  labs(title = "Comparison of European and non-European samples", x = "Scaled log likelihood ratios", 
       y = "Density", fill = "Ethnicity (based on PCA)") +
  facet_grid(rows = vars(isOfEuropeanDescent), scales = "free_y")

dev.off()

ggplot(completeCorrectedTable, aes(y = VALUE.CORRECTED, x = PGS, fill = isSampleOfInterest)) +
  geom_point(shape=20) +
  scale_fill_manual(values = c("TRUE" = "#009444", "FALSE" = "#D12F00")) +
  geom_smooth(data = completeTable, aes(x=PGS, y=VALUE.CORRECTED), method=lm, se=FALSE, color="Black", 
              size = (1 * ggplot2::.pt * 72.27/96)) +
  theme(axis.text=element_blank(),
        legend.position = "none") +
  facet_grid(rows = vars(sampleOfInterest), cols = vars(trait))