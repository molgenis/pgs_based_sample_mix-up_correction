#!/usr/bin/env Rscript

#############################################################
## c.a.warmerdam@umcg.nl - April 2020
## Plot debug figures of preliminary data
#############################################################

##############################
# Load libraries
##############################
library(tidyverse)
library(argparse)
library(data.table)
library(ROCR)
library(ggridges)

##############################
# Run stuff
##############################

traitGwasMappingPath <- "/groups/umcg-lld/tmp01/other-users/umcg-rwarmerdam/pgs_based_mixup_correction/scripts/r-scripts/pgs_based_sample_mix-up_correction/trait-gwas-mapping.txt"
basePathWithPolygenicScores <- "/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/output/PRScs/20200811/"
phenotypesFilePath <- "/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/data/lifelines/processed/pgs.phenotypes.ugli.dat"
#sampleCouplingFilePath <- "/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/data/lifelines/processed/pgs.sample-coupling-file.ugli.20201014.5120samples.txt"
sampleCouplingFilePath <- "/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/data/lifelines/processed/pgs.sample-coupling-file.ugli.20201014.perm_5120samples_51mixUps.txt"

# Calculate the scaled residuals for every combination
heightResiduals <- calculate.scaledResiduals(
  estimate = completeTable$PGS, 
  actual = completeTable$VALUE, 
  covariates = completeTable[c("AGE", "SEX")],
  responseDataType = responseDataType)

rownames(heightResiduals) <- completeTable$pheno
colnames(heightResiduals) <- completeTable$geno

meanPerRow <- rowMeans(heightResiduals)
sdPerRow <- apply(heightResiduals, 1, sd)
meanPerCol <- colMeans(heightResiduals)
sdPerCol <- apply(heightResiduals, 2, sd)

hist(meanPerRow)
hist(meanPerCol)
hist(sdPerRow)
hist(sdPerCol)

residualsScaledPerRow <- t(apply(heightResiduals, 1, function(row) {
  row <- row ^ 2
  return((row - mean(row)) / sd(row))
}))

residualsScaledPerRowDataFrame <- 
  as.data.frame.table(residualsScaledPerRow, responseName = "zScores") %>%
  inner_join(link[,c("pheno", "geno")], by = c("Var1" = "pheno"))

residualsScaledPerRowDataFrame$group <- "alternative"
residualsScaledPerRowDataFrame$group[residualsScaledPerRowDataFrame$geno == residualsScaledPerRowDataFrame$Var2] <- "null"
residualsScaledPerRowDataFrame$group <- factor(residualsScaledPerRowDataFrame$group, levels = c("null", "alternative"))

hist(residualsScaledPerRowDataFrame$zScores[residualsScaledPerRowDataFrame$group == "null"])
hist(residualsScaledPerRowDataFrame$zScores[residualsScaledPerRowDataFrame$group == "alternative"])


ggplot(residualsScaledPerRowDataFrame, aes(x=zScores, stat(density), fill=group)) +
  geom_histogram(bins = 72, alpha=.5, position="identity") +
  xlab("Scaled residuals") + ggtitle(paste0("Scaled residuals for trait '", trait, "'"))
