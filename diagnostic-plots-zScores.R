#!/usr/bin/env Rscript

#############################################################
## c.a.warmerdam@umcg.nl - October 2020
## Diagnostic plots
#############################################################

##############################
# Load libraries
##############################
library(tidyverse)
library(argparse)
library(data.table)

##############################
# Define argument parser
##############################

parser <- ArgumentParser(description='')
parser$add_argument('--trait-gwas-mapping',
                    help='path to a tab-delimited file that maps traits to the GWAS summary statistics')
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

# Load table containing paths for the plink output 
# and corresponding phenotype labels.
traitGwasMappingPath <- args$trait_gwas_mapping
traitDescriptionsTable <- fread(
  traitGwasMappingPath, 
  quote="", header=T, sep = "\t",
  col.names=c("trait", "traitDataType", "summaryStatistics", 
              "sampleSizeOfGwas", "numberOfCategories"), 
  stringsAsFactors=F)

# Get the output path
out <- args$out

# Load the phenotypes 
phenotypesFilePath <- args$phenotypes_file
phenotypesTable <- fread(phenotypesFilePath, header=T, quote="", sep="\t",
                         col.names = c("ID", "AGE", "SEX", "VALUE", "TRAIT")) %>%
  mutate(SEX = factor(SEX, levels = c("Female", "Male"))) %>%
  group_by(ID) %>%
  filter(!any(AGE < 18)) %>%
  ungroup()

# Get the link path
if (!is.null(args$sample_coupling_file)) {
  sampleCouplingFilePath <- args$sample_coupling_file
  link <- fread(sampleCouplingFilePath, stringsAsFactors=F, header=T)
}

zScoreMatrix <- matrix(nrow = nrow(link), 
                       ncol = nrow(traitDescriptionsTable), 
                       dimnames = list(link$pheno, traitDescriptionsTable$trait), NA)

# Loop trough traits
for (traitIndex in 1:nrow(traitDescriptionsTable)) {
  
  trait <- traitDescriptionsTable$trait[traitIndex]
  responseDataType <- traitDescriptionsTable$traitDataType[traitIndex]
  scaledResidualsMatrixPath <- file.path(out, trait, "scaledResidualMatrix.tsv")
  
  # Give status update
  message(paste0(traitIndex, " / ", nrow(traitDescriptionsTable), 
                 ": '", trait, "' (", responseDataType, ")."))
  
  if (!file.exists(scaledResidualsMatrixPath)){
    next
  }
  
  phenotypeTable <- phenotypesTable %>%
    filter(TRAIT == trait) %>%
    rename(pheno = ID) %>%
    inner_join(link, by="pheno")
  
  # Write scaled residuals matrix
  scaledResidualsMatrix <- fread(scaledResidualsMatrixPath) %>% 
    as.matrix(rownames = 1)
  
  scaledResidualsDataFrame <- 
    as.data.frame.table(scaledResidualsMatrix, responseName = "residuals") %>%
    inner_join(link, by = c("Var1" = "pheno")) %>%
    mutate(group = case_when(Var2 == geno ~ "null",
                             TRUE ~ "alternative"))
  
  scaledResidualsDataFrame$group <- factor(scaledResidualsDataFrame$group, levels = c("null", "alternative"))
  
  rm(scaledResidualsMatrix)
  gc()
  
  # pdf(file.path(out, trait, "unscaledResidualsHistogram.pdf"), 
  #     width=8, height = 6, useDingbats = FALSE)
  # 
  # par(xpd = NA)
  # 
  # plotResiduals(scaledResidualsDataFrame, phenotypeTable, responseDataType)
  # 
  # dev.off()
  
  residualsOfProvidedSamples <- scaledResidualsDataFrame %>% filter(group == "null")
  
  rm(scaledResidualsDataFrame)
  gc()
  
  residualsOfProvidedSamples$zScores <- scale(residualsOfProvidedSamples$residuals)
  
  zScoreMatrix[residualsOfProvidedSamples$Var1, trait] <- residualsOfProvidedSamples$zScores
}

zScoreMatrix <- zScoreMatrix[,colSums(is.na(zScoreMatrix)) < nrow(zScoreMatrix)]

zScoreDataFrame <- 
  as.data.frame.table(zScoreMatrix, responseName = "zScores") %>%
  inner_join(link, by = c("Var1" = "pheno")) %>%
  mutate(inducedMixUp = geno != original)

correctSamples <- zScoreDataFrame %>%
  filter(geno == original) %>%
  select(Var1) %>%
  distinct() %>% pull(Var1)

inducedMixUps <- zScoreDataFrame %>%
  filter(geno != original) %>%
  select(Var1) %>%
  distinct() %>% pull(Var1)

pdf(file.path(out,"residualZScoresPerTrait.correctSamples.pdf"), 
    width=8, height = 6, useDingbats = FALSE)

par(xpd = NA)

for (pheno in sample(correctSamples, size = 51)) {
    print(ggplot(zScoreDataFrame %>% filter(Var1 == pheno), 
                 aes(x=Var2, y=zScores)) +
            coord_cartesian(ylim = c(-3, 3)) +
            geom_bar(stat="identity", colour = "blue") + 
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
}

dev.off()

pdf(file.path(out,"residualZScoresPerTrait.inducedMixUps.pdf"), 
    width=8, height = 6, useDingbats = FALSE)

par(xpd = NA)

for (pheno in sample(inducedMixUps, size = 51)) {
  print(ggplot(zScoreDataFrame %>% filter(Var1 == pheno), 
        aes(x=Var2, y=zScores)) +
          coord_cartesian(ylim = c(-3, 3)) +
          geom_bar(stat="identity", colour = "red") + 
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
}

dev.off()

pdf(file.path(out,"residualZScoresPerTrait.averageOfAbsZScores.pdf"), 
    width=8, height = 6, useDingbats = FALSE)

par(xpd = NA)

zScoreSummaryTable <- zScoreDataFrame %>%
  group_by(Var2, inducedMixUp) %>%
  summarise(averageAbsoluteZScores = mean(abs(zScores), na.rm = T)) %>%
  ungroup()

print(ggplot(zScoreSummaryTable, 
             aes(x=Var2, y=averageAbsoluteZScores, fill=inducedMixUp)) +
        geom_bar(position="dodge", stat="identity") + 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))


dev.off()

meanCorrectSamples <- mean(abs(zScoreDataFrame$zScores[zScoreDataFrame$geno == zScoreDataFrame$original]), na.rm = T)
sdCorrectSamples <- sd(zScoreDataFrame$zScores[zScoreDataFrame$geno == zScoreDataFrame$original], na.rm = T)
message(paste0("Mean of abs(zScores) for correct samples: ", meanCorrectSamples))
message(paste0("SD of abs(zScores) for correct samples: ", sdCorrectSamples))

meanMixUps <- mean(abs(zScoreDataFrame$zScores[zScoreDataFrame$geno != zScoreDataFrame$original]), na.rm = T)
sdMixUps <- sd(abs(zScoreDataFrame$zScores[zScoreDataFrame$geno != zScoreDataFrame$original]), na.rm = T)
message(paste0("Mean of abs(zScores) for induced mix-ups: ", meanMixUps))
message(paste0("SD of abs(zScores) for induced mix-ups: ", sdMixUps))