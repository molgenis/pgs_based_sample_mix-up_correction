#!/usr/bin/env Rscript

#############################################################
## c.a.warmerdam@umcg.nl - August 2020
## Merges plink polygenic scores split by chromosomes.
#############################################################

##############################
# Load librairies
##############################
library(tidyverse)
library(data.table)
library(argparse)

##############################
# Define argument parser
##############################

parser <- ArgumentParser(description='')
parser$add_argument('--plink-profiles', nargs='+', required = T,
                    help='path to plink --score output files')
parser$add_argument('--out', required = T,
                    help='output path of summed profiles')
parser$add_argument('--scoresum-column', required = T,
                    help='columnname containing summable polygenic scores')

##############################
# Define functions
##############################

read.plinkProfiles <- function(plinkProfile, scoresumColumn) {
  rename_list <- list(IID = "IID", SCORESUM = scoresumColumn)
  return(fread(plinkProfile, header=T)  %>% select_(.dots = rename_list))
}

##############################
# Run
##############################

# Get the arguments
args <- parser$parse_args(commandArgs(trailingOnly = TRUE))

profileDataFrameList <- lapply(args$plink_profiles, read.plinkProfiles, args$scoresum_column)

print(sapply(profileDataFrameList, nrow))

mergedProfileDataFrame <- profileDataFrameList %>% reduce(inner_join, by = "IID")

if (length(unique(c(sapply(profileDataFrameList, nrow), mergedProfileDataFrame))) == 1) {
  stop("Profiles do not have exactly the same number of samples! Exiting...")
}

profiles <- data.frame(IID = mergedProfileDataFrame[,1], 
                       SCORESUM = apply(mergedProfileDataFrame[,-1], 1, sum))

write.table(profiles, args$out, row.names=F, col.names=T, quote=F, sep="\t")