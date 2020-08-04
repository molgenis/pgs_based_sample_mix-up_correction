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
parser$add_argument('--plink-profiles', nargs='+',
                    help='path to plink --score outputs')
parser$add_argument('--out',
                    help='output path of summed profiles')

##############################
# Define functions
##############################

read.plinkProfiles <- function(plinkProfile) {
  return(fread(plinkProfile, header=T) %>% select(IID, SCORESUM))
}

##############################
# Run
##############################

# Get the arguments
args <- parser$parse_args(commandArgs(trailingOnly = TRUE))

profileDataFrameList <- lapply(args$plink_profiles, read.plinkProfiles)

print(sapply(profileDataFrameList, nrow))

mergedProfileDataFrame <- profileDataFrameList %>% reduce(inner_join, by = "IID")

print(nrow(mergedProfileDataFrame))
print(head(mergedProfileDataFrame))

profiles <- data.frame(IID = mergedProfileDataFrame[1], 
                         SCORESUM = apply(mergedProfileDataFrame[-1], 1, sum))

write.table(profiles, args$out, row.names=F, col.names=T, quote=F, sep="\t")