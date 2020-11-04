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
parser$add_argument('--sample-coupling-file', required = FALSE,
                    help=paste0('file containing genotype sample ids in the first column',
                                'and phenotype sample ids in the second column'))
parser$add_argument('--out',
                    help='path to output directory')

##############################
# Define functions
##############################

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

# Get the link path
if (!is.null(args$sample_coupling_file)) {
  sampleCouplingFilePath <- args$sample_coupling_file
  link <- fread(sampleCouplingFilePath, stringsAsFactors=F, header=T) %>%
    filter(pheno %in% unique(phenotypesTable$ID))
}

llrMatrixOverTraits <- matrix(nrow = nrow(link), 
                    ncol = nrow(traitDescriptionsTable), 
                    dimnames = list(link$pheno, traitDescriptionsTable$trait), NA)

llrDataFrameOverTraits <- merge(data.frame(trait = traitDescriptionsTable$trait,
                                           average = NA,
                                           sd = NA),
                                data.frame(group = c("providedNotMixedUpSamples",
                                                     "providedMixedUpSamples",
                                                     "offDiagonalSamples")))

# Loop trough traits
for (traitIndex in 1:nrow(traitDescriptionsTable)) {
  
  trait <- traitDescriptionsTable$trait[traitIndex]
  responseDataType <- traitDescriptionsTable$traitDataType[traitIndex]
  llrMatrixPath <- file.path(out, trait, "logLikelihoodRatios.tsv")
  
  # Give status update
  message(paste0(traitIndex, " / ", nrow(traitDescriptionsTable), 
                 ": '", trait, "' (", responseDataType, ")."))
  
  if (!file.exists(llrMatrixPath)){
    next
  }
  
  # Write scaled residuals matrix
  llrMatrix <- fread(llrMatrixPath) %>% 
    as.matrix(rownames = 1)
  
  llrDataFrame <- 
    as.data.frame.table(llrMatrix, responseName = "logLikelihoodRatios") %>%
    inner_join(link, by = c("Var1" = "pheno")) %>%
    mutate(group = case_when(Var2 == geno ~ "null",
                             TRUE ~ "alternative"),
           inducedMixUp = original == geno)
  
  llrDataFrame$group <- factor(llrDataFrame$group, levels = c("null", "alternative"))
  
  rm(llrMatrix)
  gc()
  
  # ggplot(scaledResidualsDataFrame, aes(x=residuals, stat(density), fill=group)) +
  #   geom_histogram(bins = 72, alpha=.5, position="identity") +
  #   geom_rug(data = scaledResidualsDataFrame %>% filter(geno != original & group == "null"), aes(x=residuals), inherit.aes=F) +
  #   xlab("Unscaled residuals") + ggtitle(paste0("Unscaled residuals for trait '", trait, "'"))
  # 
  # ggsave(file.path(out, trait, "unscaledResidualsHistogram.png"), width=8, height=7)
  
  llrOfProvidedSamples <- llrDataFrame %>% filter(group == "null")

  # rm(llrDataFrame)
  # gc()
  
  llrMatrixOverTraits[llrOfProvidedSamples$Var1, trait] <- llrOfProvidedSamples$logLikelihoodRatios

  llrDataFrameOverTraits$average[llrDataFrameOverTraits$trait == trait &
                                 llrDataFrameOverTraits$group == "providedNotMixedUpSamples"] <- 
    mean(llrDataFrame[llrDataFrame$group == "null" & !llrDataFrame$inducedMixUp, "logLikelihoodRatios"], na.rm = T)
  
  llrDataFrameOverTraits$average[llrDataFrameOverTraits$trait == trait &
                                   llrDataFrameOverTraits$group == "providedMixedUpSamples"] <- 
    mean(llrDataFrame[llrDataFrame$group == "null" & llrDataFrame$inducedMixUp, "logLikelihoodRatios"], na.rm = T)
  
  llrDataFrameOverTraits$average[llrDataFrameOverTraits$trait == trait &
                                   llrDataFrameOverTraits$group == "offDiagonalSamples"] <- 
    mean(llrDataFrame[llrDataFrame$group == "alternative", "logLikelihoodRatios"], na.rm = T)
  
  llrDataFrameOverTraits$sd[llrDataFrameOverTraits$trait == trait &
                                   llrDataFrameOverTraits$group == "providedNotMixedUpSamples"] <- 
    sd(llrDataFrame[llrDataFrame$group == "null" & !llrDataFrame$inducedMixUp, "logLikelihoodRatios"], na.rm = T)
  
  llrDataFrameOverTraits$sd[llrDataFrameOverTraits$trait == trait &
                                   llrDataFrameOverTraits$group == "providedMixedUpSamples"] <- 
    sd(llrDataFrame[llrDataFrame$group == "null" & llrDataFrame$inducedMixUp, "logLikelihoodRatios"], na.rm = T)
  
  llrDataFrameOverTraits$sd[llrDataFrameOverTraits$trait == trait &
                                   llrDataFrameOverTraits$group == "offDiagonalSamples"] <- 
    sd(llrDataFrame[llrDataFrame$group == "alternative", "logLikelihoodRatios"], na.rm = T)
}

llrMatrixOverTraits <- llrMatrixOverTraits[
  ,colSums(is.na(llrMatrixOverTraits)) < nrow(llrMatrixOverTraits)]

llrDataFrameOverTraits <- 
  as.data.frame.table(llrMatrixOverTraits, responseName = "logLikelihoodRatios") %>%
  inner_join(link, by = c("Var1" = "pheno")) %>%
  mutate(inducedMixUp = geno != original)

correctSamples <- llrDataFrameOverTraits %>%
  filter(geno == original) %>%
  select(Var1) %>%
  distinct() %>% pull(Var1)

inducedMixUps <- llrDataFrameOverTraits %>%
  filter(geno != original) %>%
  select(Var1) %>%
  distinct() %>% pull(Var1)

pdf(file.path(out,"llrPerTrait.correctSamples.pdf"), 
    width=8, height = 6, useDingbats = FALSE)

par(xpd = NA)

for (pheno in sample(correctSamples, size = 51)) {
  print(ggplot(llrDataFrameOverTraits %>% filter(Var1 == pheno), 
               aes(x=Var2, y=logLikelihoodRatios)) +
          coord_cartesian(ylim = c(-3, 3)) +
          geom_bar(stat="identity", colour = "blue") + 
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
}

dev.off()

pdf(file.path(out,"llrPerTrait.inducedMixUps.pdf"), 
    width=8, height = 6, useDingbats = FALSE)

par(xpd = NA)

for (pheno in sample(inducedMixUps, size = 51)) {
  print(ggplot(llrDataFrameOverTraits %>% filter(Var1 == pheno), 
               aes(x=Var2, y=logLikelihoodRatios)) +
          coord_cartesian(ylim = c(-3, 3)) +
          geom_bar(stat="identity", colour = "red") + 
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
}

dev.off()

pdf(file.path(out,"llrPerTrait.averages.revised.pdf"), 
    width=8, height = 6, useDingbats = FALSE)

par(xpd = NA)

# llrSummaryTable <- llrDataFrameOverTraits %>%
#   group_by(Var2, inducedMixUp) %>%
#   summarise(averageLlr = mean(logLikelihoodRatios, na.rm = T)) %>%
#   ungroup()

print(ggplot(llrDataFrameOverTraits, 
             aes(x=trait, y=average, fill=group)) +
        geom_bar(position="dodge", stat="identity") + 
        geom_errorbar(aes(ymin=average-sd, ymax=average+sd), width=.2,
                      position=position_dodge(.9), alpha = 0.1) +
        coord_cartesian(ylim = c(-0.2, 0.2)) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))


dev.off()
