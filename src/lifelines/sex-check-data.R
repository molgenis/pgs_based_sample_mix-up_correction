#!/usr/bin/env Rscript

#############################################################
## c.a.warmerdam@umcg.nl - January 2021
## Plot predictive power for identifying mix-ups
#############################################################

##############################
# Define function
##############################

# Function for load sex check data for Lifelines UGLI dataset
getSexCheckData <- function(reportedSexTable) {
  ugliQc <- read_tsv("/groups/umcg-lifelines/tmp01/releases/gsa_genotypes/v1/Logs/ugli_qc_release1_samples.csv")
  gsaLinkage <- read_tsv("/groups/umcg-lifelines/tmp01/releases/gsa_linkage_files/v1/gsa_linkage_file.dat")
  
  sexCheckData <- ugliQc %>%
    inner_join(gsaLinkage, by = c("Sample_ID" = "genotyping_name")) %>%
    select(UGLI_ID, Inferred_sex) %>%
    inner_join(reportedSexTable, by = c("UGLI_ID" = "ID"))
  
  return(sexCheckData)
}

getEurSamples <- function() {
  ugliQc <- read_tsv("/groups/umcg-lifelines/tmp01/releases/gsa_genotypes/v1/Logs/ugli_qc_release1_samples.csv")
  gsaLinkage <- read_tsv("/groups/umcg-lifelines/tmp01/releases/gsa_linkage_files/v1/gsa_linkage_file.dat")
  
  eurSamples <- ugliQc %>%
    inner_join(gsaLinkage, by = c("Sample_ID" = "genotyping_name")) %>%
    filter(PCA_european) %>%
    pull(UGLI_ID)
  return(eurSamples)
}