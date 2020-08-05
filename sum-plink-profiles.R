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
# args <- parser$parse_args(c(
#   "--plink-profiles", 
#   "/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/output/PRScs/20200803//HbA1c_METAL_European.processed/1.UGLI.pgs.profile",
#   "/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/output/PRScs/20200803//HbA1c_METAL_European.processed/2.UGLI.pgs.profile",
#   "/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/output/PRScs/20200803//HbA1c_METAL_European.processed/3.UGLI.pgs.profile",
#   "/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/output/PRScs/20200803//HbA1c_METAL_European.processed/4.UGLI.pgs.profile",
#   "/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/output/PRScs/20200803//HbA1c_METAL_European.processed/5.UGLI.pgs.profile",
#   "/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/output/PRScs/20200803//HbA1c_METAL_European.processed/6.UGLI.pgs.profile",
#   "/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/output/PRScs/20200803//HbA1c_METAL_European.processed/7.UGLI.pgs.profile",
#   "/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/output/PRScs/20200803//HbA1c_METAL_European.processed/8.UGLI.pgs.profile",
#   "/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/output/PRScs/20200803//HbA1c_METAL_European.processed/9.UGLI.pgs.profile",
#   "/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/output/PRScs/20200803//HbA1c_METAL_European.processed/10.UGLI.pgs.profile",
#   "/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/output/PRScs/20200803//HbA1c_METAL_European.processed/11.UGLI.pgs.profile",
#   "/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/output/PRScs/20200803//HbA1c_METAL_European.processed/12.UGLI.pgs.profile",
#   "/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/output/PRScs/20200803//HbA1c_METAL_European.processed/13.UGLI.pgs.profile",
#   "/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/output/PRScs/20200803//HbA1c_METAL_European.processed/14.UGLI.pgs.profile",
#   "/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/output/PRScs/20200803//HbA1c_METAL_European.processed/15.UGLI.pgs.profile",
#   "/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/output/PRScs/20200803//HbA1c_METAL_European.processed/16.UGLI.pgs.profile",
#   "/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/output/PRScs/20200803//HbA1c_METAL_European.processed/17.UGLI.pgs.profile",
#   "/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/output/PRScs/20200803//HbA1c_METAL_European.processed/18.UGLI.pgs.profile",
#   "/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/output/PRScs/20200803//HbA1c_METAL_European.processed/19.UGLI.pgs.profile",
#   "/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/output/PRScs/20200803//HbA1c_METAL_European.processed/20.UGLI.pgs.profile",
#   "/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/output/PRScs/20200803//HbA1c_METAL_European.processed/21.UGLI.pgs.profile",
#   "/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/output/PRScs/20200803//HbA1c_METAL_European.processed/22.UGLI.pgs.profile",
#   "--out", 
#   "/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/output/PRScs/20200803//HbA1c_METAL_European.processed/full.UGLI.pgs.profile"
#   ))

profileDataFrameList <- lapply(args$plink_profiles, read.plinkProfiles)

print(sapply(profileDataFrameList, nrow))

mergedProfileDataFrame <- profileDataFrameList %>% reduce(inner_join, by = "IID")

if (length(unique(c(sapply(profileDataFrameList, nrow), mergedProfileDataFrame))) == 1) {
  stop("Profiles do not have exactly the same number of samples! Exiting...")
}

print(head(mergedProfileDataFrame))
print(str(mergedProfileDataFrame))

profiles <- data.frame(IID = mergedProfileDataFrame[1], 
                       SCORESUM = apply(mergedProfileDataFrame[,-1], 1, sum))

write.table(profiles, args$out, row.names=F, col.names=T, quote=F, sep="\t")