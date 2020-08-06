#!/usr/bin/env Rscript

#############################################################
## c.a.warmerdam@umcg.nl - April 2020
## Processes Lifelines phenotypes.
## Using available Lifelines phenotypes, 
## a table of processed / derived phenotypes for sample 
## mix-up correction are written with UGLI_ID identifiers.
#############################################################

##############################
# Load libraries
##############################
library(tidyverse)
library(data.table)
library(argparse)

##############################
# Get the path to the default phenotype source map
##############################

initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)
default.phenotype_sources_map <- file.path(script.basename, "phenotype-sources-map.txt")
if (length(default.phenotype_sources_map) == 0) {
  default.phenotype_sources_map = NULL
}

##############################
# Define argument parser
##############################

parser <- ArgumentParser(description='')
parser$add_argument('--phenotypes-base-path',
                    help='path to directory where tab-separated Lifelines phenotypes are stored')

parser$add_argument('--phenotype-source-map', required = FALSE, 
                    default = default.phenotype_sources_map,
                    help=paste0('tab-separated file with header and colomns indicating the phenotype name,',
                                'the column identifier and the source file.'))

parser$add_argument('--gsa-linkage-file',
                    help=paste0('file containing PSEUDOIDEXT sample ids in the first column',
                                'and UGLI_ID sample ids in the second column'))

parser$add_argument('--out',
                    help='path to a file where the derived phenotypes are to be written to.')

##############################
# Define functions
##############################

getLatestValueFromRawPhenotypeTable <- function(phenotypeSources, phenotypeTables, name) {
  columnIdentifier <- phenotypeSources$ColumnIdentifier[phenotypeSources$Name == name]
  
  return(phenotypeTables[[phenotypeSources$FileName[phenotypeSources$Name == name]]] %>%
           filter(!is.na(get(columnIdentifier))) %>%
           
           # Per PSEUDOIDEXT, get the row with the 'latest' ENCOUNTERCODE that is not NA
           group_by(PSEUDOIDEXT) %>%
           slice_max(as.numeric(ENCOUNTERCODE)) %>%
           select(PSEUDOIDEXT, !!columnIdentifier))
}

# Function for calculating eGFR
estimateGFR <- function(serum_creatine, age, female, black) {
  # Get the K constant.
  K <- if (female) 61.9 else 79.6
  # Get alpha constant.
  alpha <- if (female) -0.329 else -0.411
  
  # Calculate eGFR with base function.
  eGFR <- 141 * 
    (min(serum_creatine / K, 1) ^ alpha) * 
    (max(serum_creatine / K, 1) ^ -1.209) * 
    (0.993 ^ age)
  
  # Multiply by 1.018 if the individual is a female.
  if (female) {
    eGFR <- eGFR * 1.018
  }
  
  # Multiply by 1.159 if the individual is black.
  if (black) {
    eGFR <- eGFR * 1.159
  }
  return(eGFR)
}

##############################
# Generate phenotypes
##############################

args <- parser$parse_args(commandArgs(trailingOnly = TRUE))

args <- parser$parse_args(c(
  "--phenotypes-base-path", "/groups/umcg-lifelines/tmp01/projects/phenotypes/tab_seperated_labels/",
  "--gsa-linkage-file", "/groups/umcg-lifelines/tmp01/releases/gsa_linkage_files/v1/gsa_linkage_file.dat",
  "--out", "/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/data/lifelines/processed/UGLI.pgs.phenotypes.dat",
  "--phenotype-source-map", "/home/umcg-rwarmerdam/pgs_based_mixup_correction/scripts/r-scripts/pgs_based_sample_mix-up_correction/phenotype-source-map.txt"))

# Load GSA linkage file
gsaLinkageFile <- fread(args$gsa_linkage_file, sep="\t", header=T)

# Collect required columns
phenotypeSources <- fread(args$phenotype_source_map, sep="\t", header=T)

# For every unique file, load it, and get the respective columns.
phenotypeTables <- apply(unique(phenotypeSources[,c("Path", "FileName")]), 2,
       function(row) {
         phenotypeDataFrame <- fread(file.path(row$Path, row$FileName))
         if ("ENCOUNTERCODE" %in% colnames(phenotypeDataFrame)) {
           phenotypeDataFrame$ENCOUNTERCODE <- factor(phenotypeDataFrame$ENCOUNTERCODE, levels = c("Baseline assessment (1A)", "Second assessment (2A)"))
         }
         return(phenotypeDataFrame)
       }, USE.NAMES = TRUE, simplify = F)

# Height
height <- getLatestValueFromRawPhenotypeTable(phenotypeSources, phenotypeTables, "Height")

# BMI
BMI <- getLatestValueFromRawPhenotypeTable(phenotypeSources, phenotypeTables, "BMI")

# Diastolic blood pressure
# DBP <- getLatestValueFromRawPhenotypeTable(phenotypeSources, phenotypeTables, "Diastolic blood pressure")

# Systolic blood pressure

# Estimated GFR
estimatedGFR <- apply(phenotypes, 1, function(row) {
  if (is.na(row["geslacht"]) | (row["geslacht"] != "Female" & row["geslacht"] != "Male")) {
    return(NA)
  }
  
  return(estimateGFR(as.numeric(row["bkr"]), as.numeric(row["age_bl1"]), row["geslacht"] == "Female", FALSE))
})

# Basophil concentration

# Eosinophil concentration

# HbA1c concentration

# LDL-cholesterol

# Log-transformed triglyceride concentration

# Square root of HDL-cholesterol

# Total cholesterol concentration

# Hematocrit concentration

# Hemoglobin concentration

# Lymphocyte concentration

# Monocyte concentration

# Erythrocyte concentration

# Neutrophil concentration

# Thrombocyte concentration

# Type-2 diabetes

# Type-1 diabetes

# Asthma

# Hair colour (level of blonde (black - blonde))

# Hair colour (red vs. a level of blonde (black - blonde))

# Coronary artery disease

# Schizophrenia