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
library(MASS)
library(tidyverse)
library(argparse)
library(data.table)
library(Matrix)
library(hrbrthemes)
library(SimMultiCorrData)

##############################
# Define argument parser
##############################

parser <- ArgumentParser(description='')
parser$add_argument('--trait-gwas-mapping', required = T,
                    help='path to a tab-delimited file that maps traits to the polygenic score files.')

group <- parser$add_mutually_exclusive_group(required = T)
group$add_argument('--base-pgs-path',
                   help="path to a directory containing polygenic scores. Folder structure should be '<base-pgs-path>/<name-of-gwas-summary-statistic>/<pgs-file-name>'")
group$add_argument('--pgs-file',
                   help="path to a tab-delimited file holding all polygenic score data. Columns should represent sample IDs and individual PGSs respectively.")

parser$add_argument('--pgs-file-name',
                    help="name of files that hold polygenic scores. ignored when a combined table of polygenic scores is supplied (--pgs-file).", 
                    default = "full.UGLI.pgs.profile")

parser$add_argument('--phenotypes-file', required = T,
                    help='path to a tab-delimited file holding all processed phenotype data.')
parser$add_argument('--sample-coupling-file', required = FALSE,
                    help=paste('file containing genotype sample ids in the first column',
                               'and phenotype sample ids in the second column'))
parser$add_argument('--out',
                    help='path to output directory')

##############################
# Define constants
##############################

old <- theme_set(theme_classic())
theme_update(line = element_line(
  colour = "black", size = (0.5 / (ggplot2::.pt * 72.27/96)),
  linetype = 1, lineend = "butt", arrow = F, inherit.blank = T),
  strip.background = element_rect(colour = NA, fill = NA),
  axis.line = element_line(
    colour = "#595A5C", size = (0.5 / (ggplot2::.pt * 72.27/96)),
    linetype = 1, lineend = "butt", arrow = F, inherit.blank = T),
  axis.ticks = element_line(
    colour = "#595A5C", size = (0.5 / (ggplot2::.pt * 72.27/96)),
    linetype = 1, lineend = "butt", arrow = F, inherit.blank = T)
)

##############################
# Define functions
##############################

# Function for loading polygenic scores
loadPolygenicScores <- function(basePathWithPolygenicScores = NULL, pgsFileName = NULL, pgsMergedFile = NULL) {
  
  if (!is.null(basePathWithPolygenicScores)) {
    traitDescriptionsTable$polygenicScoreFilePath <- file.path(basePathWithPolygenicScores, 
                                                               traitDescriptionsTable$summaryStatistics, 
                                                               pgsFileName)
    
    return(bind_rows(mapply(
      function(polygenicScoreFilePath, trait) {
        message(paste0("    Loading polygenic scores from '", polygenicScoreFilePath, "'..."))
        
        # Read the PLINK polygenic score table.
        polygenicScores <- fread(
          polygenicScoreFilePath,
          header=T) %>%
          mutate(TRAIT = trait) %>%
          rename(PGS = SCORESUM, ID = IID) %>%
          dplyr::select(ID, PGS, TRAIT)
        
        return(polygenicScores)
      },
      traitDescriptionsTable$polygenicScoreFilePath, traitDescriptionsTable$trait, SIMPLIFY = F, USE.NAMES = F)))
    
  } else if (!is.null(pgsMergedFile)) {
    return(fread(pgsMergedFile, header=T, quote="", sep="\t",
                 col.names = c("ID", "PGS", "TRAIT")))
  } else {
    stop("Either 'basePathWithPolygenicScores' or 'pgsMergedFile' must not be NULL.")
  }
}

# Function to perform some type of clustering
nnit <- function(mat, clsize = 10, method=c('random','maxd', 'mind')){
  clsize.rle = rle(as.numeric(cut(1:nrow(mat), ceiling(nrow(mat)/clsize))))
  clsize = clsize.rle$lengths
  lab = rep(NA, nrow(mat))
  dmat = as.matrix(dist(mat))
  cpt = 1
  while(sum(is.na(lab)) > 0){
    lab.ii = which(is.na(lab))
    dmat.m = dmat[lab.ii,lab.ii]
    if(method[1]=='random'){
      ii = sample.int(nrow(dmat.m),1)
    } else if(method[1]=='maxd'){
      ii = which.max(rowSums(dmat.m))
    } else if(method[1]=='mind'){
      ii = which.min(rowSums(dmat.m))
    } else {
      stop('unknown method')
    }
    lab.m = rep(NA, length(lab.ii))
    lab.m[head(order(dmat.m[ii,]), clsize[cpt])] = cpt
    lab[lab.ii] = lab.m
    cpt = cpt + 1
  }
  if(any(is.na(lab))){
    lab[which(is.na(lab))] = cpt
  }
  lab
}


# Remove all individuals with red hair in the other hair phenotype
fixHairColourMutual <- function(phenotypeDataFrame, blondenessLabel, redLabel) {

  phenotypeDataFrame[phenotypeDataFrame[redLabel] == 2, blondenessLabel] <- NA_real_
  return(phenotypeDataFrame)
}

copyCorrOverCols <- function(corr) {
  copiedOverCols <- cbind(corr, corr)
  
  colnames(copiedOverCols) <- paste0(colnames(corr), 
                                     c(rep("_SIM1", times = ncol(corr)), 
                                       rep("_SIM2", times = ncol(corr))))
  return(copiedOverCols)
}

# Double matrix
enlargeCorrelationMatrix <- function(corr) {
  
  corrLarge <- cbind(corr, matrix(0, nrow = nrow(corr), ncol = ncol(corr)))
  corrLarge <- rbind(corrLarge, matrix(0, nrow = nrow(corr), ncol = ncol(corrLarge)))
  corrLarge[(nrow(corr)+1):nrow(corrLarge), (nrow(corr)+1):nrow(corrLarge)] <- corr
  
  colnames(corrLarge) <- paste0(colnames(corr), c(rep("_SIM1", times = ncol(corr)), 
                                                  rep("_SIM2", times = ncol(corr))))
  rownames(corrLarge) <- paste0(rownames(corr), c(rep("_SIM1", times = nrow(corr)), 
                                                  rep("_SIM2", times = nrow(corr))))
  
  return(corrLarge)
}

plotCorrelationMatrix <- function(corr, outPath = NULL) {
  
  pivotted <- as_tibble(corr, rownames = "rows") %>%
    pivot_longer(-rows, names_to="cols", values_to="rho")
  
  pdf(outPath, width = 8, height = 8)
  par(xpd = NA)
  
  print(ggplot(pivotted, aes(rows, cols, fill=rho)) + 
    geom_tile() +
    scale_fill_distiller(palette = "Spectral", limits = c(-1, 1)) +
    coord_equal() +
    theme(legend.position="none", axis.text.y = element_text(size = 6),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 6)))
  
  dev.off()
}

# Get mean and standard deviation for phenotypes of every trait
getPhenotypeDistributions <- function(distributionsValue, traits) {
  return(distributionsValue %>%
    filter(TRAIT %in% traits) %>%
    group_by(TRAIT) %>%
    arrange(option) %>%
    filter(1 == n()) %>% dplyr::select(TRAIT, valueMean, valueSd) %>% ungroup())
}

# Get mean and standard deviation for age covariate
getAgeDistribution <- function(distributionCov) {
  return(distributionsCov %>%
    filter(COV == "COV_AGE") %>%
    rename(TRAIT = COV) %>%
    dplyr::select(TRAIT, valueMean, valueSd) %>% ungroup())
}

# Get frequencies for the sex covariate
getSexFrequencies <- function(distributionCov) {
  return(distributionsCov %>%
    filter(COV == "COV_SEX") %>%
    rename(TRAIT = COV) %>%
    group_by(TRAIT) %>%
    arrange(option) %>%
    filter(1 != n()) %>%
    summarise(freqs = list(valueFreq[1:n()-1])) %>% ungroup())
}

# Get frequencies for the categorical covariates
getPhenotypeFrequencies <- function(distributionsValue, traits) {
  return(distributionsValue %>%
    filter(TRAIT %in% traits) %>%
    group_by(TRAIT) %>%
    arrange(option) %>%
    filter(1 != n()) %>%
    summarise(freqs = list(valueFreq[1:n()-1])) %>% ungroup())
}

getPgsDistribution <- function(distributionsPgs, pgses) {
  return(distributionsPgs %>% 
    filter(TRAIT %in% pgses) %>%
    dplyr::select(TRAIT, c("valueMean" = "pgsMean", 
                           "valueSd" = "pgsSd")) %>% ungroup())
}

hyperbolicFactorTan <- function(rho, x) {
  return(sqrt(tanh(atanh(rho^2) * x)))
}

correlationComparison <- function(corrA, corrB, outPath = NULL) {
  plotCorrelationMatrix(corrA, outPath = paste0(outPath, "_targetCorrelations.pdf"))
  plotCorrelationMatrix(corrB, outPath = paste0(outPath, "_actualCorrelations.pdf"))
  
  bound <- inner_join(
    as_tibble(corrA, rownames = "rows") %>%
      pivot_longer(-rows, names_to="cols", values_to="targetRho"),
    as_tibble(corrB, rownames = "rows") %>%
      pivot_longer(-rows, names_to="cols", values_to="actualRho"), 
    by = c("rows", "cols")) %>%
    mutate(type = case_when(
      rows == cols ~ "diagonal",
      startsWith(rows, "VALUE") & startsWith(cols, "VALUE") ~ "pheno - pheno",
      (startsWith(rows, "PGS") & startsWith(cols, "VALUE")) 
      | (startsWith(rows, "PGS") & startsWith(cols, "VALUE")) ~ "pgs - pheno",
      startsWith(rows, "PGS") & startsWith(cols, "PGS") ~ "pgs - pgs",
      startsWith(rows, "COV") | startsWith(cols, "COV") ~ "cov - pgs, cov - pheno"))
  
  pdf(paste0(outPath, "_scatterComparison.pdf"), width = 8, height = 8)
  par(xpd = NA)
  
  print(ggplot(bound, aes(targetRho, actualRho, colour = type)) +
    geom_point() +
    coord_equal())
  
  dev.off()
}

# Generate new phenotypes
phenotypeSimulation <- function(baseCorr, 
                                distributionsPgs, distributionsCov, distributionsPheno, 
                                clustering,
                                outDir, nTraits, explainedVarianceFactor) {
  
  # Make output folder
  folder <- paste0(outDir, "_", format(Sys.Date(), "%Y%m%d"))
  dir.create(file.path(folder))
  
  # Adapt correlation matrix
  pgsSel <- startsWith(rownames(baseCorr), "PGS_")
  valueSel <- startsWith(colnames(baseCorr), "VALUE_")
  covsSel <- rownames(baseCorr) == c("COV_SEX", "COV_AGE")
  
  # Get all variables that are in the matrix in each categories
  pgses <- rownames(baseCorr)[pgsSel]
  traits <- rownames(baseCorr)[valueSel]
  covs <- rownames(baseCorr)[covsSel]
  
  # Get distributions and frequencies
  distValueTibble <- getPhenotypeDistributions(distributionsPheno, traits)
  freqsValueTibble <- getPhenotypeFrequencies(distributionsPheno, traits)
  distAgeTibble <- getAgeDistribution(distributionsCov)
  freqsSexTibble <- getSexFrequencies(distributionCov)
  distPgsTibble <- getPgsDistribution(distributionsPgs, pgses)
  
  # We will need to assign probabilities to each of the phenotypes according to their AUC
  
  nTraitsToSample <- as.numeric(nTraits) %% nrow(clustering)
  
  nDups <- as.numeric(nTraits) %/% nrow(clustering)
  
  # Adapted matrix
  longBaseCorr <- as_tibble(baseCorr, rownames = "rows") %>%
    pivot_longer(-rows, names_to="cols", values_to="rho") %>%
    separate(rows, c("rowType", "rowName"), sep = "_") %>%
    separate(cols, c("colType", "colName"), sep = "_") %>%
    left_join(clustering %>% rename(rowClust = traitClust), by = c("rowName" = "traitName")) %>%
    left_join(clustering %>% rename(colClust = traitClust), by = c("colName" = "traitName"))
    
  # Copy the base correlations x times
  withDuplicates <- mapply(function(copyIndex) {
    return(longBaseCorr %>% 
      mutate(copyIndex = as.integer(copyIndex)))
  }, c(1:(nDups+1)), SIMPLIFY = F)
  
  # From the last interation, remove n rows if necessary
  withDuplicates[[nDups+1]] <- withDuplicates[[nDups+1]] %>%
    filter(is.na(rowClust) | (rowClust * 5) <= nTraitsToSample, 
           is.na(colClust) | (colClust * 5) <= nTraitsToSample)
  
  # Merge the remaining parts.
  completeCorrLong <- bind_rows(withDuplicates) %>%
    mutate(
      newRowNames = case_when(
        rowType == "PGS" | rowType == "VALUE" ~ paste0(rowType, "_", rowName, "_SIM", copyIndex),
        TRUE ~ paste0(rowType, "_", rowName)),
      newColNames = case_when(
        colType == "PGS" | colType == "VALUE" ~ paste0(colType, "_", colName, "_SIM", copyIndex),
        TRUE ~ paste0(colType, "_", colName))) %>%
    dplyr::select(-copyIndex) %>%
    distinct()
  
  # Inflate correlations between PGSs and phenotypes by a certain amount
  inflatedCorrLong <- completeCorrLong %>%
    mutate(rho = case_when(
      colType == "VALUE" & rowType == "PGS" ~ hyperbolicFactorTan(rho, explainedVarianceFactor),
      rowType == "VALUE" & colType == "PGS" ~ hyperbolicFactorTan(rho, explainedVarianceFactor),
      TRUE ~ rho))
  
  adaptedMatrixLong <- inflatedCorrLong %>% 
    dplyr::select(newRowNames, newColNames, rho) %>%
    pivot_wider(id_cols = newRowNames, 
                names_from = newColNames, values_from = rho,
                values_fill = 0)
  
  corrLongSelection <- completeCorrLong %>% 
    filter(rowType %in% c("PGS", "VALUE")) %>%
    dplyr::select(rowName, rowType, newRowNames) %>% distinct() %>%
    mutate(combinedRowName = paste0(rowType, "_", rowName)) 
  
  distValueTibble2 <- distValueTibble %>%
    inner_join(corrLongSelection, by = c("TRAIT" = "combinedRowName")) %>% 
    dplyr::select(-TRAIT) %>% rename(TRAIT = newRowNames)
  distPgsTibble2 <- distPgsTibble %>%
    inner_join(corrLongSelection, by = c("TRAIT" = "combinedRowName")) %>% 
    dplyr::select(-TRAIT) %>% rename(TRAIT = newRowNames)
  freqsValueTibble2 <- freqsValueTibble %>%
    inner_join(corrLongSelection, by = c("TRAIT" = "combinedRowName")) %>% 
    dplyr::select(-TRAIT) %>% rename(TRAIT = newRowNames)
  
  # Bind the distributions and the frequencies
  distTibble <- bind_rows(distAgeTibble, distValueTibble2, distPgsTibble2)
  freqsTibble <- bind_rows(freqsSexTibble, freqsValueTibble2)
  orderedTraits <- c(freqsTibble$TRAIT, distTibble$TRAIT)
  
  adaptedMatrix <- adaptedMatrixLong %>% dplyr::select(-newRowNames) %>% as.matrix()
  rownames(adaptedMatrix) <- adaptedMatrixLong$newRowNames
  
  adaptedPdMatrix <- nearPD(adaptedMatrix, corr = T, base.matrix = T)$mat
  rownames(adaptedPdMatrix) <- rownames(adaptedMatrix)
  colnames(adaptedPdMatrix) <- colnames(adaptedMatrix)
  
  n_cat <- nrow(freqsTibble)
  n_cont <- nrow(distTibble)
  
  corr <- adaptedPdMatrix[orderedTraits, orderedTraits]
  
  valid_out <- valid_corr(k_cat = n_cat, k_cont = n_cont, 
                          method = "Polynomial", 
                          means = distTibble$valueMean, 
                          vars = distTibble$valueSd^2,
                          skews = rep(0, n_cont),
                          skurts = rep(0, n_cont),
                          fifths = rep(0, n_cont),
                          sixths = rep(0, n_cont),
                          marginal = freqsTibble$freqs, rho = corr, n = 10000, seed = 1234)
  
  corr_out <- rcorrvar(k_cat = n_cat, k_cont = n_cont, 
                       method = "Polynomial", 
                       means = distTibble$valueMean, 
                       vars = distTibble$valueSd^2,
                       skews = rep(0, n_cont),
                       skurts = rep(0, n_cont),
                       fifths = rep(0, n_cont),
                       sixths = rep(0, n_cont),
                       marginal = freqsTibble$freqs, rho = corr, n = 50000, seed = 1234)
  
  ordinalNew <- corr_out[["ordinal_variables"]]
  contNew <- corr_out[["continuous_variables"]]
  colnames(ordinalNew) <- freqsTibble$TRAIT
  colnames(contNew) <- distTibble$TRAIT
  
  newPhenotypeDataFrame <- cbind(ordinalNew, contNew)

  hcTable <- freqsTibble %>% 
    filter(startsWith(TRAIT, "VALUE_")) %>% 
    separate(TRAIT, c("type", "name", "sim"), sep = "_") %>%
    mutate(TRAIT = paste0("VALUE_", name, "_", sim)) %>%
    filter(name %in% c("Blondeness of hair", "Red hair colour")) %>% 
    dplyr::select(name, sim, TRAIT) %>% pivot_wider(names_from = name, values_from = TRAIT) %>%
    arrange(sim)
    
  for (i in 1:nrow(hcTable)) {
    newPhenotypeDataFrame <- fixHairColourMutual(
      newPhenotypeDataFrame, 
      hcTable$`Blondeness of hair`[i], 
      hcTable$`Red hair colour`[i])
  }

  correlationsActual <- cor(newPhenotypeDataFrame, use = "pairwise.complete.obs")

  correlationComparison(adaptedMatrix, correlationsActual, outPath = file.path(folder, "correlationComparison"))
  
  newPhenotypeDataFrame$ID <- rownames(newPhenotypeDataFrame)
  
  # We should make a new description table
  # We should remove the prefixes from the names
  
  # We need to inner_join the phenotypes table with the rownames table to double rows from the phenotypes table
  # When this is required.
  phenotypesTable <- newPhenotypeDataFrame %>%
    dplyr::select(ID, starts_with("COV_"), starts_with("VALUE_")) %>%
    rename(SEX = COV_SEX, AGE = COV_AGE) %>%
    pivot_longer(cols = starts_with("VALUE_"), names_to = "TRAIT", values_to = "VALUE", names_prefix = "VALUE_")
  
  polygenicScoresTable <- newPhenotypeDataFrame %>%
    dplyr::select(ID, starts_with("PGS_")) %>%
    pivot_longer(cols = starts_with("PGS_"), names_to = "TRAIT", values_to = "PGS", names_prefix = "PGS_")
  
  traitDescriptionsTable <- phenotypesTable %>%
    dplyr::select(-ID, -SEX, -AGE) %>%
    rename(trait = TRAIT) %>%
    group_by(trait) %>%
    summarise(
      traitDataType = case_when(
        all(paste0("VALUE_", trait) %in% freqsTibble$TRAIT) ~ "ordinal",
        TRUE ~ "continuous"),
      numberOfCategories = case_when(
        all(traitDataType == "ordinal") ~ length(na.omit(unique(VALUE))),
        TRUE ~ 1L),
      traitDataType = case_when(all(numberOfCategories == 2) ~ "binary", 
                                all(numberOfCategories > 2) ~ "ordinal", 
                                TRUE ~ traitDataType)
    )
  
  write.table(corr, file.path(folder, "targetCorrelations.txt"), sep = "\t", row.names = T, col.names = T, quote = F)
  write.table(correlationsActual, file.path(folder, "actualCorrelations.txt"), sep = "\t", row.names = T, col.names = T, quote = F)
  write.table(phenotypesTable, file.path(folder, "phenotypesTable.txt"), sep = "\t", row.names = F, col.names = T, quote = F)
  write.table(polygenicScoresTable, file.path(folder, "polygenicScoresTable.txt"), sep = "\t", row.names = F, col.names = T, quote = F)
  write.table(traitDescriptionsTable, file.path(folder, "traitDescriptionTable.txt"), sep = "\t", row.names = F, col.names = T, quote = F)
}

##############################
# Run
##############################

args <- parser$parse_args(c(
  "--trait-gwas-mapping",
  "/groups/umcg-lld/tmp01/other-users/umcg-rwarmerdam/pgs_based_mixup_correction/scripts/r-scripts/pgs_based_sample_mix-up_correction/data/lifelines/trait-gwas-mapping.txt",
  "--base-pgs-path",
  "/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/output/PRScs/20201120/",
  "--phenotypes-file", 
  "/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/data/lifelines/processed/pgs.phenotypes_20201215.ugli.dat",
  "--out", "."
  ))

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
polygenicScoresTable <- loadPolygenicScores(
  basePathWithPolygenicScores = args$base_pgs_path, pgsFileName = args$pgs_file_name, 
  pgsMergedFile = args$pgs_file)

# Get the output path
out <- args$out

message(strwrap(prefix = " ", initial = "", paste(
  "Loading phenotype table:\n", args$phenotypes_file)))

# Load the phenotypes 
phenotypesFilePath <- args$phenotypes_file
phenotypesTable <- fread(phenotypesFilePath, header=T, quote="", sep="\t") %>%
  distinct() %>%
  rename_all(recode, "UGLI_ID" = "ID") %>%
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

# Link the phenotypes table to the sample coupling file
phenotypesTable <- phenotypesTable %>%
  rename(pheno = ID) %>%
  inner_join(link, by="pheno")

# Merge the PGSs with actual phenotype data
completeLongTable <- phenotypesTable %>%
  inner_join(polygenicScoresTable, by = c("geno" = "ID", "TRAIT" = "TRAIT")) %>%
  rename(ID = pheno)

completeTable <- completeLongTable %>% as_tibble() %>%
  select(ID, VALUE, TRAIT, PGS, AGE, SEX) %>%
  filter(!is.na(VALUE) & !is.na(PGS))

covsTable <- completeTable %>% 
  select(ID, SEX, AGE) %>%
  group_by(ID) %>%
  summarise(SEX = head(SEX, n=1), AGE = head(AGE, n=1)) %>%
  mutate(SEX = as.integer(SEX)) %>%
  ungroup()

distributionsValue <- completeTable %>% inner_join(traitDescriptionsTable %>% select(trait, traitDataType), by = c("TRAIT" = "trait")) %>%
  select(ID, VALUE, TRAIT, PGS, traitDataType) %>%
  distinct() %>%
  ungroup() %>%
  group_by(TRAIT) %>%
  mutate(option = case_when(traitDataType != "continuous" ~ VALUE),
         totalN = n()) %>%
  ungroup() %>%
  group_by(TRAIToption = paste(TRAIT, as.character(option), sep = '_')) %>%
  summarise(valueMean = mean(VALUE, na.rm = T),
            valueSd = sd(VALUE, na.rm = T),
            valueN = n(),
            valueFreq = unique(n() / totalN),
            TRAIT = paste0("VALUE_", head(TRAIT, 1)),
            option = head(as.character(option), 1))

distributionsCov <- covsTable %>%
  pivot_longer(c(SEX, AGE), names_to = "COV", values_to = "VALUE") %>%
  group_by(COV) %>%
  mutate(option = case_when(COV == "SEX" ~ VALUE,
                            COV == "AGE" ~ NA_integer_),
         totalN = n()) %>%
  ungroup() %>%
  group_by(TRAIToption = paste(COV, as.character(option), sep = '_')) %>%
  summarise(valueMean = mean(VALUE, na.rm = T),
            valueSd = sd(VALUE, na.rm = T),
            valueN = n(),
            valueFreq = unique(n() / totalN),
            COV = paste0("COV_", head(COV, 1)),
            option = head(as.character(option), 1))

distributionsPgs <- completeTable %>% inner_join(traitDescriptionsTable %>% select(trait, traitDataType), by = c("TRAIT" = "trait")) %>%
  select(ID, VALUE, TRAIT, PGS) %>%
  distinct() %>%
  ungroup() %>%
  group_by(TRAIT) %>%
  mutate(totalN = n()) %>%
  summarise(pgsMean = mean(PGS, na.rm = T),
            pgsSd = sd(PGS, na.rm = T),
            valueN = n(),
            valueFreq = unique(n() / totalN),
            TRAIT = paste0("PGS_", head(TRAIT, 1)))

write.table(distributionsValue, "descVal_20210622.txt", sep = "\t", row.names = F, col.names = T, quote = F)
write.table(distributionsCov, "descCov_20210622.txt", sep = "\t", row.names = F, col.names = T, quote = F)
write.table(distributionsPgs, "descPgs_20210622.txt", sep = "\t", row.names = F, col.names = T, quote = F)

phenoWide <- completeTable %>% 
  pivot_wider(c(ID), values_from = c(VALUE, PGS), names_from = TRAIT) %>%
  inner_join(covsTable, by = c("ID")) %>%
  rename(COV_AGE = AGE, COV_SEX = SEX)
  
fullMatrix <- phenoWide %>% select(-ID) %>% as.matrix()
rownames(fullMatrix) <- phenoWide %>% pull(ID)

correlationMatrix <- cor(fullMatrix, use = "pairwise.complete.obs")
write.table(correlationMatrix, "cormat_20210622.txt", sep = "\t", row.names = T, col.names = T, quote = F)

correlationMatrix <- as.matrix(read.table(
  "~/Documents/projects/pgs_based_mixup_correction-ugli/jobs/pgs-simulations/cormat_20210622.txt", 
  sep = "\t", quote = "", header = T, check.names = F))
distributionsPgs <- read.table(
  "~/Documents/projects/pgs_based_mixup_correction-ugli/jobs/pgs-simulations/descPgs_20210622.txt", 
  sep = "\t", quote = "", header = T)
distributionsCov <- read.table(
  "~/Documents/projects/pgs_based_mixup_correction-ugli/jobs/pgs-simulations/descCov_20210622.txt", 
  sep = "\t", quote = "", header = T)
distributionsPheno <- read.table(
  "~/Documents/projects/pgs_based_mixup_correction-ugli/jobs/pgs-simulations/descVal_20210622.txt", 
  sep = "\t", quote = "", header = T)

toRemove <- c("VALUE_Type 1 diabetes", "PGS_Type 1 diabetes")
correlationMatrix <- correlationMatrix[
  !(rownames(correlationMatrix) %in% toRemove), 
  !(colnames(correlationMatrix) %in% toRemove)]

plotCorrelationMatrix(corr = correlationMatrix, 
                      out = file.path(out, "baseCorrelations.pdf"))

# The correlation between Red hair colour and Blondeness of hair is NA, since there is no overlap
# Here we assume red hair colour to not be correlated with your blondeness of hair; Getting 
# red hair depends on genetics for red hair colour, irrespective of how blonde your hair would be.
correlationMatrix[is.na(correlationMatrix)] <- 0

# If we want to perform a fine step-based approach of 5 or 10 traits, we should make an order of the
# phenotypes to determine how they should be added.

# Perhaps we should cluster the phenotypes based on correlation matrix and the AUC that the cluster
# has.
# We should list the correlations between each 

corrLong <- as_tibble(correlationMatrix, rownames = "rows") %>%
  pivot_longer(-rows, names_to="cols", values_to="rho") %>%
  separate(rows, c("rowType", "rowName"), sep = "_") %>%
  separate(cols, c("colType", "colName"), sep = "_") %>%
  filter(colType == "VALUE", rowType == "PGS")

mat <- t(as.matrix(corrLong %>% 
                     pivot_wider(names_from = "colName", values_from = "rho") %>% 
                     dplyr::select(-rowName, -rowType, -colType)))

clsize <- 5

lab = nnit(mat, clsize, method="maxd")

corrLong2 <- corrLong %>% 
  left_join(tibble(colName = rownames(mat), colClust = lab), by = "colName") %>%
  left_join(tibble(rowName = rownames(mat), rowClust = lab), by = "rowName") %>%
  mutate(colName = fct_reorder(colName, colClust),
         rowName = fct_reorder(rowName, rowClust))
  # pivot_wider(names_from = "colName", values_from = "rho") %>% 
  
main <- ggplot(corrLong2, aes(x = colName, y = rowName, fill = rho)) + 
  geom_tile() +
  scale_fill_distiller(palette = "Spectral", limits = c(-1, 1)) +
  coord_equal() +
  theme_minimal() +
  theme(legend.position="none", axis.text.y = element_text(size = 6),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 6))

main

pc <- ggplot(corrLong2, aes(x = colName, y=0, fill = colClust)) + geom_tile() + 
  scale_y_discrete(position="right") +
  theme_minimal() +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) +
  xlab(NULL) + ylab(NULL)

clustering <- tibble(traitName = rownames(mat), traitClust = lab)
# 
# %>%
#   filter(colType == "VALUE", rowType == "PGS") %>%
#   group_by(colName) %>%
#   filter(rowName == colName) %>%
#   group_by() %>%
#   mutate(meanRho = mean(rho), sdRho = sd(rho)) %>%
#   ungroup() %>%
#   mutate(prob = dnorm(rho, mean = meanRho, sd = sdRho),
#          distance = )

sweepFactor <- expand_grid(explainedVarianceFactor = c(0.5), 
                           nTraits = c(10, 25, 50, 75, 100))

result <- apply(sweepFactor, 1, function(row) {
  outDir <- file.path(out, paste0("simulationDataset_nTraits", row["nTraits"], 
         "_evf", row["explainedVarianceFactor"]))
  
  phenotypeSimulation(
    baseCorr = correlationMatrix,
    distributionsPgs = distributionsPgs, 
    distributionsCov = distributionsCov, 
    distributionsPheno = distributionsPheno, 
    outDir = outDir, 
    clustering = clustering,
    nTraits = row["nTraits"],
    explainedVarianceFactor = row["explainedVarianceFactor"])
  
  return(outDir)
})
