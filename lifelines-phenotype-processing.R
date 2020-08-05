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
# Define argument parser
##############################

parser <- ArgumentParser(description='')
parser$add_argument('--phenotypes_base_path',
                    help='path to directory where tab-separated Lifelines phenotypes are stored')

parser$add_argument('--phenotype_source_map', required = FALSE, 
                    default = paste0(dirname(sys.frame(1)$ofile), "/phenotype-sources-map.txt"),
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



##############################
# Generate phenotypes
##############################

args <- parser$parse_args(commandArgs(trailingOnly = TRUE))

print(args)
# args <- parser$parse_args(c("--phenotypes_base_path", ""))