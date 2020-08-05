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

getHeight <- function() {
  
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

# Collect required columns

# Height


# BMI

# Diastolic blood pressure

# Systolic blood pressure

# Estimated GFR

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

# Reticulocyte concentration

# Neutrophil concentration

# Thrombocyte concentration

# Type-2 diabetes

# Type-1 diabetes

# Asthma

# Hair colour (level of blonde (black - blonde))

# Hair colour (red vs. a level of blonde (black - blonde))

# Coronary artery disease

# Schizophrenia