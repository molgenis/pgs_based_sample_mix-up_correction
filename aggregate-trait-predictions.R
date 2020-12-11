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
library(pROC)

##############################
# Define argument parser
##############################

parser <- ArgumentParser(description='')
# parser$add_argument('--debug', action='store_true', 
#                     dest="debug", help="write intermediate llr and residuals files.")
# parser$add_argument('--base-fit-model-path',
#                     help='path to a directory to write fitted model parmaters to, or to load fitted models from.')
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
# Run
##############################
# args <- parser$parse_args(c("--trait-gwas-mapping", "/groups/umcg-lld/tmp01/other-users/umcg-rwarmerdam/pgs_based_mixup_correction/scripts/r-scripts/pgs_based_sample_mix-up_correction/trait-gwas-mapping.txt",
#                             "--sample-coupling-file", "/home/umcg-rwarmerdam/pgs_based_mixup_correction-ugli/data/lifelines/processed/pgs.sample-coupling-file.ugli.20201120.10080samples.txt",
#                             "--base-pgs-path", "/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/output/PRScs/20201120/",
#                             "--phenotypes-file", "/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/data/lifelines/processed/pgs.phenotypes.ugli.dat",
#                             "--out", "/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/output/sample-swap-prediction/20200811.test/",
#                             "--llr-bayes-method", "NA", "50",
#                             "--base-fit-model-path", "/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/output/sample-swap-prediction/fitted-paramters/"))
args <- parser$parse_args(commandArgs(trailingOnly = TRUE))

message(strwrap(prefix = " ", initial = "", paste(
  "Loading trait-gwas-mapping:\n", args$trait_gwas_mapping)))

# Load table containing paths for the plink output 
# and corresponding phenotype labels.
traitGwasMappingPath <- args$trait_gwas_mapping
traitDescriptionsTable <- fread(
  traitGwasMappingPath, 
  quote="", header=T, sep = "\t",
  col.names=c("trait", "traitDataType", "summaryStatistics", 
              "sampleSizeOfGwas", "numberOfCategories"), 
  stringsAsFactors=F)

message(strwrap(prefix = " ", initial = "", paste(
  "Loading polygenic scores from:\n", args$trait_gwas_mapping)))

# Get the paths to the polygenic scores.
basePathWithPolygenicScores <- args$base_pgs_path
traitDescriptionsTable$polygenicScoreFilePath <- file.path(basePathWithPolygenicScores, 
                                                 traitDescriptionsTable$summaryStatistics, 
                                                 "full.UGLI.pgs.profile")

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
                  "Not using special sample coupling table"))
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

traitDescriptionsTable <- traitDescriptionsTable %>%
  mutate(confinedAuc = NA_real_,
         matrixWideAuc = NA_real_,
         pValue = NA_real_,
         traitOutputDir = NA_character_,
         modelBasePath = modelBasePath)
  
# Loop trough traits

for (traitIndex in 1:nrow(traitDescriptionsTable)) {
  
  trait <- traitDescriptionsTable$trait[traitIndex]

  #traitFileName <- paste(traitIndex, gsub(" ", "_", trait), sep = ".")
  traitFileName <- trait
  traitDescriptionsTable[traitIndex, "traitOutputDir"] <- paste(traitIndex, gsub(" ", "_", trait), sep = ".")

  logLikelihoodRatios <- as.matrix(fread(file.path(out, traitFileName, "/logLikelihoodRatios.tsv")), rownames = 1)
            
  likelihoodRatioDifferenceTest <- t.test(
    diag(logLikelihoodRatios), 
    logLikelihoodRatios[lower.tri(logLikelihoodRatios) | upper.tri(logLikelihoodRatios)],
    alternative = "less")
  
  print(likelihoodRatioDifferenceTest)
  
  if (likelihoodRatioDifferenceTest$p.value <= likelihoodRatioDifferenceAlpha) {
    
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
  
  llrDataFrame <- 
    as.data.frame.table(logLikelihoodRatios, responseName = "logLikelihoodRatios", stringsAsFactors = FALSE) %>%
    inner_join(link, by = c("Var1" = "pheno"))
  
  rm(logLikelihoodRatios)
  gc()
  
  llrDataFrame$group <- "alternative"
  llrDataFrame$group[llrDataFrame$original == llrDataFrame$Var2] <- "null"
  llrDataFrame$group <- ordered(llrDataFrame$group, levels = c("null", "alternative"))
  
  matrixWideAuc <- auc(
    llrDataFrame$group, llrDataFrame$logLikelihoodRatios)
  
  message(paste0("Calculated overall AUC: ", matrixWideAuc))
  
  permutationTestDataFrame <- llrDataFrame %>%
    filter(geno == Var2)
  
  traitDescriptionsTable[traitIndex, "matrixWideAuc"] <- as.double(matrixWideAuc)
  traitDescriptionsTable[traitIndex, "pValue"] <- likelihoodRatioDifferenceTest$p.value
  
  if (predictingInducedMixUps && "alternative" %in% permutationTestDataFrame$group) {
    
    confinedAuc <- auc(
      permutationTestDataFrame$group, 
      permutationTestDataFrame$logLikelihoodRatios)
    
    message(paste0("Confined AUC: ", confinedAuc))
    traitDescriptionsTable[traitIndex, "confinedAuc"] <- as.double(confinedAuc)
  }
  
  # Convert matrix to data frame for further processing
  llrDataFrame <- 
    as.data.frame.table(aggregatedLlrMatrix, responseName = "logLikelihoodRatios") %>%
    inner_join(link, by = c("Var1" = "pheno"))
  
  llrDataFrame$group <- "alternative"
  llrDataFrame$group[llrDataFrame$original == llrDataFrame$Var2] <- "null"
  llrDataFrame$group <- ordered(llrDataFrame$group, levels = c("null", "alternative"))
  
  newAuc <- auc(
    llrDataFrame$group, llrDataFrame$logLikelihoodRatios)
  message(paste0("Calculated overall AUC: ", newAuc))
}

# Write the result values
# write.table(traitDescriptionsTable, file.path(out, "outputStatisticsPerTrait.tsv"),
#             sep="\t", col.names = T, row.names = T, quote = F)
# 
# write.table(aggregatedLlrMatrix, file.path(out, "aggregatedLogLikelihoodRatiosMatrix.tsv"),
#             sep="\t", col.names = T, row.names = T, quote = F)
# 
# write.table(aggregatedNumberOfTraits, file.path(out, "aggregatedNumberOfTraitsMatrix.tsv"),
#             sep="\t", col.names = T, row.names = T, quote = F)

print(dim(aggregatedLlrMatrix))

# aggregatedLlrMatrix <- aggregatedLlrMatrix[apply(aggregatedNumberOfTraits > 0, 1, any), 
#                                            apply(aggregatedNumberOfTraits > 0, 2, any)]

overallOutputStatistics <- data.frame(name = "overallStatistics")

likelihoodRatioDifferenceTest <- t.test(
  diag(aggregatedLlrMatrix), 
  aggregatedLlrMatrix[lower.tri(aggregatedLlrMatrix) | upper.tri(aggregatedLlrMatrix)],
  alternative = "less")

overallOutputStatistics$pValue <- likelihoodRatioDifferenceTest$p.value

# Convert matrix to data frame for further processing
llrDataFrame <- 
  as.data.frame.table(aggregatedLlrMatrix, responseName = "logLikelihoodRatios") %>%
  inner_join(link, by = c("Var1" = "pheno"))

# write.table(llrDataFrame, file.path(out, "/aggregatedLogLikelihoodRatiosDataFrame.tsv"),
#             sep="\t", col.names = T, row.names = F, quote = F)

rm(aggregatedLlrMatrix)
gc()

llrDataFrame$group <- "alternative"
llrDataFrame$group[llrDataFrame$original == llrDataFrame$Var2] <- "null"
llrDataFrame$group <- ordered(llrDataFrame$group, levels = c("null", "alternative"))

matrixWideAuc <- auc(
  llrDataFrame$group, llrDataFrame$logLikelihoodRatios)

message(paste0("Calculated overall AUC: ", matrixWideAuc))
overallOutputStatistics$matrixWideAuc <- as.double(matrixWideAuc)

# pdf(file.path(out, "ROCcurve_matrixWide.pdf"))
# par(xpd = NA)

# roc(
#   llrDataFrame$group ~ llrDataFrame$logLikelihoodRatios, plot=TRUE,
#   print.auc=TRUE,col="green",lwd =4,legacy.axes=TRUE,main="ROC Curves")
# 
# dev.off()

# Confine ourselves to the diagonal
permutationTestDataFrame <- llrDataFrame %>%
  filter(geno == Var2)

if ("alternative" %in% permutationTestDataFrame$group) {
  confinedAuc <- auc(
    permutationTestDataFrame$group, 
    permutationTestDataFrame$logLikelihoodRatios)
  
  message(paste0("Confined AUC: ", confinedAuc))
  overallOutputStatistics$confinedAuc <- as.double(confinedAuc)
  
  # pdf(file.path(out, "ROCcurve_diagonal.pdf"))
  # par(xpd = NA)
  # 
  # roc(
  #   permutationTestDataFrame$group ~ permutationTestDataFrame$logLikelihoodRatios, 
  #   plot=TRUE, print.auc=TRUE,col="green",lwd =4,legacy.axes=TRUE,main="ROC Curves")
  # 
  # dev.off()
}

# write.table(overallOutputStatistics, file.path(out, "overallOutputStatistics.tsv"),
#             sep="\t", col.names = T, row.names = T, quote = F)
# 
# ggplot(llrDataFrame, aes(x=logLikelihoodRatios, stat(density), fill=group)) +
#   geom_histogram(bins = 32, alpha=.5, position="identity") +
#   xlab("Log likelihood ratios") + ggtitle(paste0("LR overall"))
# 
# ggsave(file.path(out, "/likelihoodRatioHistogram.png"), width=8, height=7)
# 
# ggplot(llrDataFrame, aes(x=group, y=logLikelihoodRatios)) +
#   geom_boxplot() + ggtitle(paste0("Log likelihood ratio distributions comparison overall"))
# 
# ggsave(file.path(out, "/likelihoodRatioBoxplot.png"), width=8, height=7)
# 
# numberOfTraits <- 
#   as.data.frame.table(aggregatedNumberOfTraits, responseName = "numberOfTraits")
# 
# phenoSamples <- unique(llrDataFrame[llrDataFrame$geno != llrDataFrame$original, "Var1"])

# phenoSamples <- sample(unique(lrProducts[lrProducts$geno == lrProducts$original, "Var1"]), size = 26)
# 
# sampledLrProducts <- lrProducts %>% 
#   #inner_join(numberOfTraits, by=c("Var1", "Var2")) %>%
#   mutate(colourGroup = case_when(
#     Var2 == geno ~ 1,
#     Var2 == original ~ 2,
#     TRUE ~ 3)) %>%
#   filter(Var1 %in% phenoSamples)
# 
# cols = c("2" = "blue", "1" = "red", "3" = "black")
# 
# ggplot(sampledLrProducts %>% filter(is.finite(logLikelihoodRatios)), 
#        aes(y = Var1)) +
#   
#   # Add ridge lines
#   geom_density_ridges(
#     aes(x = logLikelihoodRatios, 
#         point_alpha = as.numeric(colourGroup != 3),
#         point_color = colourGroup),
#     jittered_points = TRUE,
#     position = position_points_jitter(width = 0, height = 0),
#     point_shape = '|', point_size = 3, alpha = 0.7) +
#   scale_point_color_continuous(low = "#0072B2", high = "#D55E00") +
#   #scale_discrete_manual(values = cols, aesthetics = "point_color") +
#   #scale_discrete_manual("point_color", values = cols) +
# 
#   # Add annotation with number of traits, and the number of filtered likelihood ratios
#   # geom_text(
#   #   data=sampledLrProducts %>% group_by(Var1) %>%
#   #     summarise(numberOfTraits = median(numberOfTraits),
#   #               numberOfFiltered = sum(!is.finite(logLikelihoodRatios))),
#   #   position=position_nudge(y=0.64), colour="red", size=3.5,
#   #   hjust = "inward", x = 0,
#   #   aes(label = sprintf("n traits: %d, n filtered: %d", numberOfTraits, numberOfFiltered))) +
# 
#   # Set theme
#   theme_minimal() +
#   theme(legend.position = "none")
# 
# ggsave("~/pgs_based_mixup_correction-ugli/output/sample-swap-prediction/20200811.20201012.efi-discretization.30/ridges.png", width=8, height=20)

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

# sampleTraitsSummary <- sampledLrProducts %>% group_by(Var1) %>%
#       summarise(numberOfTraits = median(numberOfTraits),
#                 numberOfFiltered = sum(!is.finite(logLikelihoodRatios)))
# 
# lrProducts_Prev <- 
#   as.data.frame.table(aggregatedLlrMatrix, responseName = "logLikelihoodRatios") %>%
#   inner_join(link, by = c("Var1" = "pheno"))
# 
# lrProducts_Prev$group <- "alternative"
# lrProducts_Prev$group[lrProducts_Prev$original == lrProducts_Prev$Var2] <- "null"
# lrProducts_Prev$group <- factor(lrProducts_Prev$group, c("alternative", "null"))
