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
library(readr)
library(data.table)
library(argparse)

##############################
# Get the path to the default phenotype source map
##############################

initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)
default.phenotype_sources_map <- file.path(script.basename, "phenotype-source-map.txt")
if (length(default.phenotype_sources_map) == 0) {
  default.phenotype_sources_map = NULL
}

ENCOUNTERCODE_LEVELS <- c(
  "Baseline assessment (1A)", 
  "Follow-up questionnaire (1B)", 
  "Follow-up questionnaire (1C)", 
  "Second assessment (2A)")

VMID_LEVELS <- c(
  "Baseline assessment, birth questionnaire",
  "Baseline assessment, questionnaire 1",
  "Baseline assessment, questionnaire 2",
  "Baseline assessment, questionnaire 3",
  "Follow-up questionnaire (1B)",
  "Follow-up questionnaire (1C)", 
  "Second assessment, questionnaire 1",
  "Second assessment, questionnaire 2")

##############################
# Define argument parser
##############################

parser <- ArgumentParser(description='')
parser$add_argument('--phenotype-source-map', required = FALSE, 
                    default = default.phenotype_sources_map,
                    help=paste0('tab-separated file with header and colomns indicating the phenotype name,',
                                'the column identifier and the source file.'))

parser$add_argument('--gsa-linkage-file',
                    help=paste0('file containing PSEUDOIDEXT sample ids in the first column',
                                'and UGLI_ID sample ids in the second column'))

parser$add_argument('--out',
                    help='prefix where the derived phenotypes are to be written to.')

##############################
# Define functions
##############################

loadPhenotypeTables <- function(filePath) {

  phenotypeDataFrame <- fread(filePath,
    header=T, quote="", sep="\t", fill=T)
  
  message(paste0("Dimensions of '", filePath, "': ", nrow(phenotypeDataFrame), " x ", ncol(phenotypeDataFrame)))
  
  if ("ENCOUNTERCODE" %in% colnames(phenotypeDataFrame)) {
    
    phenotypeDataFrame$ENCOUNTERCODE <- factor(
      phenotypeDataFrame$ENCOUNTERCODE, 
      levels = ENCOUNTERCODE_LEVELS)
  }
  
  if ("VMID" %in% colnames(phenotypeDataFrame)) {
    
    phenotypeDataFrame$VMID <- factor(
      phenotypeDataFrame$VMID, 
      levels = VMID_LEVELS)
  }
  
  return(phenotypeDataFrame %>% mutate(PSEUDOIDEXT = as.character(PSEUDOIDEXT)))
}

getVmidFromEncountercode <- function(encountercodes) {
  encounterVmidMap <- c(
    "Baseline assessment (1A)" = "Baseline assessment, questionnaire 1",
    "Follow-up questionnaire (1B)" = "Follow-up questionnaire (1B)",
    "Follow-up questionnaire (1C)" = "Follow-up questionnaire (1C)", 
    "Second assessment (2A)" = "Second assessment, questionnaire 1"
  )
  
  return(factor(unname(encounterVmidMap[encountercodes]), levels = VMID_LEVELS))
}

getCorrectionTable <- function(phenotypeSources, phenotypeTables) {
  vmidAgeNameMap <- c(
    "Baseline assessment, birth questionnaire" = "Age 1A1",
    "Baseline assessment, questionnaire 1" = "Age 1A1",
    "Baseline assessment, questionnaire 2" = "Age 1A2",
    "Baseline assessment, questionnaire 3" = "Age 1A3",
    "Follow-up questionnaire (1B)" = "Age 1B",
    "Follow-up questionnaire (1C)" = "Age 1C", 
    "Second assessment, questionnaire 1" = "Age 2A1",
    "Second assessment, questionnaire 2" = "Age 2A2"
  )
  
  # Get columns that contain ages for every VMID
  vmidAgeColumnMap <- sapply(vmidAgeNameMap, function(name) phenotypeSources$ColumnIdentifier[phenotypeSources$Name == name], USE.NAMES=T)
  # Get the column identifier for the sex of individuals
  sexColumnIdentifier <- phenotypeSources$ColumnIdentifier[phenotypeSources$Name == "Sex"]

  correctionColumnIdentifiers <- c(vmidAgeColumnMap, sexColumnIdentifier)
  tableNamesWithCorrectionValues <- unique(phenotypeSources$filePath[phenotypeSources$Name %in% c(vmidAgeNameMap, "Sex")])
  
  if (length(tableNamesWithCorrectionValues) != 1) {
    stop("could not retrieve age and sex from a single table! Exiting")
  }
  
  correctionTable <- phenotypeTables[[tableNamesWithCorrectionValues]] %>%
    select(PSEUDOIDEXT, !!correctionColumnIdentifiers) %>%
    pivot_longer(
      cols = names(vmidAgeColumnMap),
      names_to = "VMID",
      values_to = "AGE",
      values_drop_na = TRUE
    ) %>%
    mutate(VMID = factor(VMID, 
                         levels = VMID_LEVELS)) %>%
    rename_all(recode, GESLACHT = "SEX")
  
  return(correctionTable)
}

getBlondenessOfHair <- function(phenotypeSources, phenotypeTables, correctionTable) {
  columnIdentifier <- phenotypeSources$ColumnIdentifier[phenotypeSources$Name == "Natural hair colour now"]
  return(phenotypeTables[[phenotypeSources$filePath[phenotypeSources$Name == "Natural hair colour now"]]] %>%
           mutate(
             VALUE = case_when(get(columnIdentifier) == "Fair" ~ 2L,
                               get(columnIdentifier) == "Brown" ~ 1L,
                               get(columnIdentifier) == "Black" ~ 0L,
                               TRUE ~ NA_integer_)) %>%
           filter(!is.na(VALUE)) %>%
           inner_join(correctionTable, by = c("PSEUDOIDEXT", "VMID")) %>%
           group_by(PSEUDOIDEXT) %>%
           slice(1) %>%
           select(PSEUDOIDEXT, AGE, SEX, VALUE))
}

getRedHairValues <- function(phenotypeSources, phenotypeTables, correctionTable) {
  columnIdentifier <- phenotypeSources$ColumnIdentifier[phenotypeSources$Name == "Natural hair colour now"]
  return(phenotypeTables[[phenotypeSources$filePath[phenotypeSources$Name == "Natural hair colour now"]]] %>%
           mutate(
             VALUE = case_when(
               get(columnIdentifier) %in% c("Fair", "Brown", "Black") ~ 0L,
               get(columnIdentifier) == "Red" ~ 1L,
               TRUE ~ NA_integer_)) %>%
           
           filter(!is.na(VALUE)) %>%
           inner_join(correctionTable, by = c("PSEUDOIDEXT", "VMID")) %>%
           group_by(PSEUDOIDEXT) %>%
           slice(1) %>%
           select(PSEUDOIDEXT, AGE, SEX, VALUE))
}

getSchizophreniaValues <- function(phenotypeSources, phenotypeTables, correctionTable) {
  columnIdentifier <- phenotypeSources$ColumnIdentifier[phenotypeSources$Name == "Schizophrenia"]
  return(phenotypeTables[[phenotypeSources$filePath[phenotypeSources$Name == "Schizophrenia"]]] %>%
           mutate(
             VALUE = case_when(get(columnIdentifier) == "Yes" ~ 1L,
                               TRUE ~ 0L)) %>%
           filter(!is.na(VALUE)) %>%
           inner_join(correctionTable, by = c("PSEUDOIDEXT", "VMID")) %>%
           group_by(PSEUDOIDEXT) %>%
           slice(1) %>%
           select(PSEUDOIDEXT, AGE, SEX, VALUE))
}

getCoronaryArteryDiseaseValues <- function(phenotypeSources, phenotypeTables, correctionTable) {
  
  # Get tables indicating MI, or surgery for every question
  # Get Baseline MI
  columnIdentifierMi <- phenotypeSources$ColumnIdentifier[phenotypeSources$Name == "Heart attack"]
  
  baselineMi <- phenotypeTables[[phenotypeSources$filePath[phenotypeSources$Name == "Heart attack"]]] %>%
    mutate(
      VALUE = case_when(get(columnIdentifierMi) == "Yes" ~ 1L,
                        get(columnIdentifierMi) == "No" ~ 0L,
                        TRUE ~ NA_integer_)) %>%
    filter(!is.na(VALUE)) %>%
    inner_join(correctionTable, by = c("PSEUDOIDEXT", "VMID")) %>%
    select(PSEUDOIDEXT, AGE, SEX, VALUE, VMID)
  
  # Get follow up Mi
  columnIdentifierMiFollowUp <- phenotypeSources$ColumnIdentifier[phenotypeSources$Name == "Since last questionnaire: heart attack"]
  
  followUpMi <- phenotypeTables[[phenotypeSources$filePath[phenotypeSources$Name == "Since last questionnaire: heart attack"]]] %>%
    mutate(
      VALUE = case_when(get(columnIdentifierMiFollowUp) == "Yes" ~ 1L,
                        get(columnIdentifierMiFollowUp) == "No" ~ 0L,
                        TRUE ~ NA_integer_)) %>%
    filter(!is.na(VALUE)) %>%
    inner_join(correctionTable, by = c("PSEUDOIDEXT", "VMID")) %>%
    select(PSEUDOIDEXT, AGE, SEX, VALUE, VMID)
  
  # Get baseline intervention
  columnIdentifierInterventionBaseline <- phenotypeSources$ColumnIdentifier[phenotypeSources$Name == "Balloon angioplasty/bypass surgery"]
  
  baselineIntervention <- phenotypeTables[[phenotypeSources$filePath[phenotypeSources$Name == "Balloon angioplasty/bypass surgery"]]] %>%
    mutate(
      VALUE = case_when(get(columnIdentifierInterventionBaseline) == "Yes" ~ 1L,
                        get(columnIdentifierInterventionBaseline) == "No" ~ 0L,
                        TRUE ~ NA_integer_)) %>%
    filter(!is.na(VALUE)) %>%
    inner_join(correctionTable, by = c("PSEUDOIDEXT", "VMID")) %>%
    select(PSEUDOIDEXT, AGE, SEX, VALUE, VMID)
  
  columnIdentifierInterventionFollowUp <- phenotypeSources$ColumnIdentifier[phenotypeSources$Name == "Balloon angioplasty since last questionnaire"]
  
  followUpIntervention <- phenotypeTables[[phenotypeSources$filePath[phenotypeSources$Name == "Balloon angioplasty since last questionnaire"]]] %>%
    mutate(
      VALUE = case_when(get(columnIdentifierInterventionFollowUp) == "Yes" ~ 1L,
                        TRUE ~ NA_integer_)) %>%
    filter(!is.na(VALUE)) %>%
    inner_join(correctionTable, by = c("PSEUDOIDEXT", "VMID")) %>%
    select(PSEUDOIDEXT, AGE, SEX, VALUE, VMID)
  
  cad <- bind_rows(list(
    "HEALTH19" = baselineMi, 
    "HEALTH100B1" = followUpMi, 
    "HEALTH23" = baselineIntervention, 
    "HEALTH101M1" = followUpIntervention), .id = "ID") %>%
    group_by(PSEUDOIDEXT) %>%
    slice_max(VALUE, with_ties=FALSE) %>%
    select(PSEUDOIDEXT, AGE, SEX, VALUE)
  
  return(cad)
}

getDerivedBinaryValues <- function(phenotypeSources, phenotypeTables, correctionTable, name, encountercode) {
  columnIdentifier <- phenotypeSources$ColumnIdentifier[phenotypeSources$Name == name]
  return(phenotypeTables[[phenotypeSources$filePath[phenotypeSources$Name == name]]] %>%
           mutate(
             VMID = factor(getVmidFromEncountercode(encountercode),
                                    levels = VMID_LEVELS),
             VALUE = case_when(get(columnIdentifier) == "Yes" ~ 1L,
                               get(columnIdentifier) == "No" ~ 0L,
                               TRUE ~ NA_integer_)) %>%
           filter(!is.na(VALUE)) %>%
           inner_join(correctionTable, by = c("PSEUDOIDEXT", "VMID")) %>%
           group_by(PSEUDOIDEXT) %>%
           slice(1) %>%
           select(PSEUDOIDEXT, AGE, SEX, VALUE))
}

getBloodPressureValues <- function(phenotypeSources, phenotypeTables, correctionTable, name) {
  columnIdentifier <- phenotypeSources$ColumnIdentifier[phenotypeSources$Name == name]
  return(phenotypeTables[[phenotypeSources$filePath[phenotypeSources$Name == name]]] %>%
           mutate(VMID = factor(getVmidFromEncountercode("Baseline assessment (1A)"), 
                                         levels = VMID_LEVELS),
                  VALUE = get(columnIdentifier)) %>%
           filter(!is.na(VALUE)) %>%
           inner_join(correctionTable, by = c("PSEUDOIDEXT", "VMID")) %>%
           group_by(PSEUDOIDEXT) %>%
           slice(1) %>%
           select(PSEUDOIDEXT, AGE, SEX, VALUE))
}

getBloodTraitPhenotypesTable <- function(phenotypeSources, phenotypeTables, correctionTable) {
  # Method for getting a table with all blood traits. 
  # (PGSes calculated with summary statistics from Vuckovic et al. (2020))
  # We have to grab the variables from below
  bloodTraits <- c("Basophilic Granulocytes",
                   "Eosinophil concentration",
                   "Lymphocytes",
                   "Leukocytes",
                   "Monocytes",
                   "Erythrocytes",
                   "Neutrophil Granulocytes",
                   "Thrombocytes",
                   "Hematocrit",
                   "Hemoglobin")
  
  names <- c(bloodTraits,
             "Pregnancy Baseline",
             "Pregnancy Followup",
             "Cancer, cured",
             "Cancer",
             "Cancer Followup",
             "Renal failure",
             "Liver cirrhosis")
  
  # 0b. List the names as column identifiers in a vector
  
  columnIdentifiers <- sapply(names, function(name) {
    phenotypeSources$ColumnIdentifier[phenotypeSources$Name == name]
    })
  
  # 1. First, get a list of extracted tables with the requested columns.
  
  extractedPhenotypeTableList <- extractPhenotypeTables(phenotypeSources = phenotypeSources,
                                                        phenotypeTables = phenotypeTables,
                                                        names = names)
  
  # 2. Then, get a table containing if the individual likely has cancer
  
  # Table with cancer at baseline.
  hasCancerAtBaselineTable <- extractedPhenotypeTableList[c("Cancer, cured", "Cancer")] %>% 
    reduce(inner_join, by = c("PSEUDOIDEXT", "VMID")) %>%
    filter(VMID == "Baseline assessment, questionnaire 1") %>%
    mutate(hasCancer = case_when(
      # Set to one (1) if cancer is reported AND the individual was not declared cancer free, 
      # OR if the individual reports to have had cancer, and is not declared cancer free.
      (get(columnIdentifiers["Cancer"]) == "Yes" 
      & get(columnIdentifiers["Cancer, cured"]) != "Yes")
      | get(columnIdentifiers["Cancer, cured"]) == "No" ~ 1L,
      
      # Set to zero (0) if cancer was not reported, OR the individual reports to have been
      # declared cancer free.
      get(columnIdentifiers["Cancer"]) == "No" 
      | get(columnIdentifiers["Cancer, cured"]) == "Yes" ~ 0L,
      
      # For other cases, set to NA
      TRUE ~ NA_integer_)) %>%
    filter(!is.na(hasCancer)) %>%
    select(PSEUDOIDEXT, VMID, hasCancer)
  
  # Table with cancer since last questionnaire.
  hasCancerAtFollowupTable <- extractedPhenotypeTableList[["Cancer Followup"]] %>%
    mutate(hasCancer = case_when(
      # Set to one (1) if cancer started since last questionnaire
      get(columnIdentifiers["Cancer Followup"]) == "Yes" ~ 1L,
      # Set to one (0) if cancer has not started since last questionnaire
      get(columnIdentifiers["Cancer Followup"]) == "No" ~ 0L,
      # Set to NA if neither No nor Yes was answered.
      TRUE ~ NA_integer_)) %>%
    filter(!is.na(hasCancer)) %>%
    select(PSEUDOIDEXT, VMID, hasCancer)
  
  # Bind tables.
  hasCancerTable <- bind_rows(hasCancerAtBaselineTable, hasCancerAtFollowupTable)

  # Remove unnecessary variables
  rm(hasCancerAtBaselineTable)
  rm(hasCancerAtFollowupTable)
  gc()
  
  # 3. Get a table containing pregnancy at 1A
  
  pregnancyAtBaseline <- extractedPhenotypeTableList[["Pregnancy Baseline"]] %>% 
    mutate(isPregnant = case_when(
      get(columnIdentifiers["Pregnancy Baseline"]) == "Yes" ~ 1L,
      get(columnIdentifiers["Pregnancy Baseline"]) == "No" ~ 0L,
      TRUE ~ NA_integer_)) %>%
    filter(!is.na(isPregnant)) %>%
    select(PSEUDOIDEXT, VMID, isPregnant)
  
  pregnancyAtFollowup <- extractedPhenotypeTableList[["Pregnancy Followup"]] %>% 
    mutate(isPregnant = case_when(
      get(columnIdentifiers["Pregnancy Followup"]) == "Yes" ~ 1L,
      get(columnIdentifiers["Pregnancy Followup"]) == "No" ~ 0L,
      TRUE ~ NA_integer_)) %>%
    filter(!is.na(isPregnant)) %>%
    select(PSEUDOIDEXT, VMID, isPregnant)
  
  # Bind tables.
  isPregnantTable <- bind_rows(pregnancyAtBaseline, pregnancyAtFollowup)
  
  # Remove unnecessary variables
  rm(pregnancyAtBaseline)
  rm(pregnancyAtFollowup)
  gc()
  
  # 4. Get a table containing renal failure
  
  renalFailureTable <- extractedPhenotypeTableList[["Renal failure"]] %>%
    mutate(hasRenalFailure = case_when(
      get(columnIdentifiers["Renal failure"]) == "Yes" ~ 1L,
      TRUE ~ NA_integer_)) %>%
    filter(!is.na(hasRenalFailure)) %>%
    select(PSEUDOIDEXT, hasRenalFailure)
  
  # 5. Get a table containing liver cirrhosis
  
  liverCirrhosisTable <- extractedPhenotypeTableList[["Liver cirrhosis"]] %>%
    mutate(hasLiverCirrhosis = case_when(
      get(columnIdentifiers["Liver cirrhosis"]) == "Yes" ~ 1L,
      TRUE ~ NA_integer_)) %>%
    filter(!is.na(hasLiverCirrhosis)) %>%
    select(PSEUDOIDEXT, hasLiverCirrhosis)
  
  # 4 Merge all tables (on PSUEDOIDEXT and VMID if applicable), 
  # and make sure one of the VMID (if applicable) is retained
  
  bloodTraitPhenotypesTable <- extractedPhenotypeTableList[bloodTraits] %>% 
    reduce(full_join, by = c("PSEUDOIDEXT", "VMID")) %>%
    left_join(hasCancerTable, by = c("PSEUDOIDEXT", "VMID")) %>%
    left_join(isPregnantTable, by = c("PSEUDOIDEXT", "VMID")) %>%
    left_join(renalFailureTable, by = "PSEUDOIDEXT") %>%
    left_join(liverCirrhosisTable, by = "PSEUDOIDEXT") %>%
    filter(
      # Remove rows with Platelet > 1000 * 10^(9) /L
      !(get(columnIdentifiers["Thrombocytes"]) > 1000)

      # Remove rows with WBC count > 200 * 10^(9) /L
      & !(get(columnIdentifiers["Leukocytes"]) > 200)

      # Remove rows with Hemoglobin > 20g /dL
      & !(get(columnIdentifiers["Hemoglobin"]) > 20 * 0.6206)

      # Remove rows with Hematocrit > 0.6 (L/L)
      & !(get(columnIdentifiers["Hematocrit"]) > 0.6)

      & (is.na(hasCancer) | hasCancer != 1)
      & (is.na(isPregnant) | isPregnant != 1)
      & (is.na(hasRenalFailure) | hasRenalFailure != 1)
      & (is.na(hasLiverCirrhosis) | hasLiverCirrhosis != 1)) %>%
    select(PSEUDOIDEXT, VMID, all_of(columnIdentifiers[bloodTraits])) %>%
    inner_join(correctionTable, by = c("PSEUDOIDEXT", "VMID")) %>%
    
    # Per PSEUDOIDEXT, get the row with the 'latest' VMID that is not NA
    group_by(PSEUDOIDEXT) %>%
    slice_max(as.numeric(VMID)) %>%
    select(-VMID)
  
  return(bloodTraitPhenotypesTable)
}

extractPhenotypeTables <- function(phenotypeSources, phenotypeTables, names) {
  
  phenotypeTableList <- sapply(names, function(name) {
    columnIdentifier <- phenotypeSources$ColumnIdentifier[phenotypeSources$Name == name]
    
    return(phenotypeTables[[phenotypeSources$filePath[phenotypeSources$Name == name]]] %>%
      mutate(VMID = if("VMID" %in% colnames(.)) VMID else getVmidFromEncountercode(ENCOUNTERCODE)) %>%
      filter(!is.na(get(columnIdentifier))) %>%
      select(PSEUDOIDEXT, VMID, all_of(columnIdentifier)))
    }, simplify = FALSE, USE.NAMES = TRUE)
  
  return(phenotypeTableList)
}

getLatestValueFromRawPhenotypeTable <- function(phenotypeSources, phenotypeTables, correctionTable, name, 
                                                correctionColumns = c("AGE", "SEX")) {
  
  columnIdentifier <- phenotypeSources$ColumnIdentifier[phenotypeSources$Name == name]
  
  return(phenotypeTables[[phenotypeSources$filePath[phenotypeSources$Name == name]]] %>%
           mutate(VALUE = get(columnIdentifier),
                  VMID = if("VMID" %in% colnames(.)) VMID else getVmidFromEncountercode(ENCOUNTERCODE)) %>%
           filter(!is.na(VALUE) & VALUE != "") %>%
           inner_join(correctionTable, by = c("PSEUDOIDEXT", "VMID")) %>%
           
           # Per PSEUDOIDEXT, get the row with the 'latest' VMID that is not NA
           group_by(PSEUDOIDEXT) %>%
           slice_max(as.numeric(VMID)) %>%
           select(PSEUDOIDEXT, all_of(correctionColumns), VALUE))
}

getEduYears <- function(phenotypeSources, phenotypeTables, correctionTable) {
  return(getLatestValueFromRawPhenotypeTable(
    phenotypeSources, phenotypeTables, correctionTable, "Highest completed education") %>%
      mutate(VALUE = case_when(
          VALUE == "Junior general secondary education (such as MAVO, (M)ULO, MBO-short, VMBO-t)" ~ 10L,
          VALUE == "lower or preparatory secondary vocational education (such as LTS, LEAO, LHNO, VMBO)" ~ 10L,
          VALUE == "No education (did not finish primary school)" ~ 1L,
          VALUE == "Primary education (primary school, special needs primary school)" ~ 7L,
          VALUE == "Secondary vocational education or work-based learning pathway (such as MBO-long, VWO, Atheneum, Gymnasium, INAS)" ~ 15L,
          VALUE == "Senior general secondary education, pre-university secondary education (such as HAVO, VWO, Atheneum, gymnasium, HBS, MMS" ~ 13L,
          VALUE == "Higher vocational education (such as HBO, HTS, HEAO, doctoral university education, Bachelor's)" ~ 22L,
          VALUE == "University education" ~ 22L
        )) %>%
      filter(!is.na(VALUE) & AGE >= 30))
}

getValueFromBloodTraitPhenotypesTable <- function(
  bloodTraitPhenotypesTable, name) {
  
  # Rename the value of interest to 'VALUE', and remove NA's.
  # Return only the requested columns.
  return(bloodTraitPhenotypesTable %>%
    rename(VALUE = !!name) %>%
    filter(!is.na(VALUE)) %>%
    select(PSEUDOIDEXT, AGE, SEX, VALUE))
}

getSquareRootOfHdlCholesterol <- function(phenotypeSources, phenotypeTables, correctionTable) {
  return(getLatestValueFromRawPhenotypeTable(
    phenotypeSources, phenotypeTables, correctionTable, "HDL Cholesterol") %>%
      mutate(VALUE = sqrt(VALUE)))
}

getLogTransformedTriglycerideConcentration <- function(phenotypeSources, phenotypeTables, correctionTable) {
  return(getLatestValueFromRawPhenotypeTable(
    phenotypeSources, phenotypeTables, correctionTable, "Triglycerides") %>%
      mutate(VALUE = log(VALUE)))
}

getLogTransformedEosinophilConcentration <- function(bloodTraitPhenotypesTable) {
  return(getValueFromBloodTraitPhenotypesTable(
    bloodTraitPhenotypesTable, "Eosinophil concentration") %>%
    mutate(VALUE = log10(VALUE + 0.01)))
}

getLogTransformedNeutrophilConcentration <- function(bloodTraitPhenotypesTable) {
  return(getValueFromBloodTraitPhenotypesTable(
    bloodTraitPhenotypesTable, "Neutrophil Granulocytes") %>%
      mutate(VALUE = log10(VALUE + 0.01)))
}

getLogTransformedBasophilConcentration <- function(bloodTraitPhenotypesTable) {
  return(getValueFromBloodTraitPhenotypesTable(
    bloodTraitPhenotypesTable, "Basophilic Granulocytes") %>%
      mutate(VALUE = log10(VALUE + 0.01)))
}

getLogTransformedMonocyteConcentration <- function(bloodTraitPhenotypesTable) {
  return(getValueFromBloodTraitPhenotypesTable(
    bloodTraitPhenotypesTable, "Monocytes") %>%
      mutate(VALUE = log10(VALUE + 0.01)))
}

getLogTransformedLymphocyteConcentration <- function(bloodTraitPhenotypesTable) {
  return(getValueFromBloodTraitPhenotypesTable(
    bloodTraitPhenotypesTable, "Lymphocytes") %>%
      mutate(VALUE = log10(VALUE + 0.01)))
}

getEstimatedGfr <- function(phenotypeSources, phenotypeTables, correctionTable) {
  creatinine <- getLatestValueFromRawPhenotypeTable(
    phenotypeSources, phenotypeTables, correctionTable, "Creatinine")

  return(creatinine %>%
    filter(!is.na(SEX) & SEX %in% c("Male", "Female")) %>%
    mutate(VALUE = log(estimateGfr(VALUE, AGE, SEX == "Female", FALSE))) %>%
    filter(!is.na(VALUE)))
}

# Function for calculating eGFR
estimateGfr <- function(serum_creatine, age, female, black) {
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

plotPhenotype <- function(name, tbl) {
  #frequencies <- table(tbl$VALUE, useNA="ifany")
  phenotypeHistogram <- ggplot(tbl, aes(x=VALUE)) + geom_histogram() + labs(x=name, y="Count")
}

##############################
# Generate phenotypes
##############################

args <- parser$parse_args(commandArgs(trailingOnly = TRUE))

args <- parser$parse_args(c("--gsa-linkage-file", "/groups/umcg-lifelines/tmp01/releases/gsa_linkage_files/v1/gsa_linkage_file.dat","--out", "/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/data/lifelines/processed/UGLI.pgs.phenotypes_tmp.dat","--phenotype-source-map", "/home/umcg-rwarmerdam/pgs_based_mixup_correction/scripts/r-scripts/pgs_based_sample_mix-up_correction/phenotype-source-map.txt"))

message("Started.")
# Load GSA linkage file
gsaLinkageTable <- fread(args$gsa_linkage_file, sep="\t", header=T) %>% 
  mutate(PSEUDOIDEXT = as.character(PSEUDOIDEXT))

# Collect required columns
phenotypeSources <- fread(args$phenotype_source_map, sep="\t", header=T) %>%
  mutate(filePath = file.path(Path, FileName))

message("Loading phenotype tables.")
# For every unique file, load it, and get the respective columns.
phenotypeTables <- sapply(unique(phenotypeSources$filePath),
       loadPhenotypeTables, USE.NAMES = TRUE, simplify = F)
message("Loaded all tables.")

# Getting correction table
correctionTable <- getCorrectionTable(phenotypeSources, phenotypeTables)
message("Created correction table.")

bloodTraitPhenotypesTable <- getBloodTraitPhenotypesTable(phenotypeSources = phenotypeSources,
                                                          phenotypeTables = phenotypeTables,
                                                          correctionTable = correctionTable)
message("Created blood trait phenotypes table.")

traitList <- list()

message("Processing phenotypes...")

# Height
message("    Heigth...")
traitList[["Height"]] <- getLatestValueFromRawPhenotypeTable(
  phenotypeSources, phenotypeTables, correctionTable, "Height")

# BMI
message("    BMI...")
traitList[["BMI"]] <- getLatestValueFromRawPhenotypeTable(
  phenotypeSources, phenotypeTables, correctionTable, "BMI")

# Diastolic blood pressure
message("    Diastolic blood pressure...")
traitList[["Diastolic blood pressure"]] <- getBloodPressureValues(
  phenotypeSources, phenotypeTables, correctionTable, "Diastolic blood pressure")

# Systolic blood pressure
message("    Systolic blood pressure...")
traitList[["Systolic blood pressure"]] <- getBloodPressureValues(
  phenotypeSources, phenotypeTables, correctionTable, "Systolic blood pressure")

# Estimated GFR
message("    eGFR...")
traitList[["eGFR"]] <- getEstimatedGfr(
  phenotypeSources, phenotypeTables, correctionTable)

# Basophil concentration
message("    Basophilic granulocytes...")
traitList[["Basophilic Granulocyte concentration"]] <- getLogTransformedBasophilConcentration(
  bloodTraitPhenotypesTable)

# Eosinophil concentration
message("    Eosinophil concentration...")
traitList[["Eosinophil concentration"]] <- getLogTransformedEosinophilConcentration(
  bloodTraitPhenotypesTable)

# HbA1c concentration
message("    HbA1c concentration")
traitList[["HbA1c"]] <- getLatestValueFromRawPhenotypeTable(
  phenotypeSources, phenotypeTables, correctionTable, "HbA1c")

# LDL-cholesterol
message("    LDL-cholesterol...")
traitList[["LDL cholesterol"]] <- getLatestValueFromRawPhenotypeTable(
  phenotypeSources, phenotypeTables, correctionTable, "LDL Cholesterol")

# Log-transformed triglyceride concentration
message("    Triglyceride concentration...")
traitList[["Triglyceride concentration"]] <- getLogTransformedTriglycerideConcentration(
  phenotypeSources, phenotypeTables, correctionTable)

# Square root of HDL-cholesterol
message("    HDL-cholesterol...")
traitList[["HDL cholesterol"]] <- getSquareRootOfHdlCholesterol(
  phenotypeSources, phenotypeTables, correctionTable)

# Total cholesterol concentration
message("    Total cholesterol...")
traitList[["Total cholesterol"]] <- getLatestValueFromRawPhenotypeTable(
  phenotypeSources, phenotypeTables, correctionTable, "Cholesterol")

# Hematocrit concentration
message("    Hematocrit concentration...")
traitList[["Hematocrit concentration"]] <- getValueFromBloodTraitPhenotypesTable(
  bloodTraitPhenotypesTable, "Hematocrit")

# Hemoglobin concentration
message("    Hemoglobin concentration...")
traitList[["Hemoglobin concentration"]] <- getValueFromBloodTraitPhenotypesTable(
  bloodTraitPhenotypesTable, "Hemoglobin")

# Lymphocyte concentration
message("    Lymphocyte concentration...")
traitList[["Lymphocyte concentration"]] <- getLogTransformedLymphocyteConcentration(
  bloodTraitPhenotypesTable)

# Monocyte concentration
message("    Monocyte concentration...")
traitList[["Monocyte concentration"]] <- getLogTransformedMonocyteConcentration(
  bloodTraitPhenotypesTable)

# Erythrocyte concentration
message("    Erythrocyte concentration...")
traitList[["Erythrocyte concentration"]] <- getValueFromBloodTraitPhenotypesTable(
  bloodTraitPhenotypesTable, "Erythrocytes")

# Neutrophil concentration
message("    Neutrophil concentration...")
traitList[["Neutrophil concentration"]] <- getLogTransformedNeutrophilConcentration(
  bloodTraitPhenotypesTable)

# Thrombocyte concentration
message("    Thrombocyte concentration...")
traitList[["Thrombocyte concentration"]] <- getValueFromBloodTraitPhenotypesTable(
  bloodTraitPhenotypesTable, "Thrombocytes")

# Type-2 diabetes
message("    Type 2 diabetes...")
traitList[["Type 2 diabetes"]] <- getDerivedBinaryValues(
  phenotypeSources, phenotypeTables, correctionTable, "Type 2 diabetes", "Second assessment (2A)")

# Type-1 diabetes
message("    Type 1 diabetes...")
traitList[["Type 1 diabetes"]] <- getDerivedBinaryValues(
  phenotypeSources, phenotypeTables, correctionTable, "Type 1 diabetes", "Second assessment (2A)")

# Asthma
message("    Asthma...")
traitList[["Asthma"]] <- getDerivedBinaryValues(
  phenotypeSources, phenotypeTables, correctionTable, "Asthma", "Baseline assessment (1A)")

# Hair colour (level of blonde (black - blonde))
message("    Blondeness of hair...")
traitList[["Blondeness of hair"]] <- getBlondenessOfHair(
  phenotypeSources, phenotypeTables, correctionTable)

# Hair colour (red vs. a level of blonde (black - blonde))
message("    Red hair colour...")
traitList[["Red hair colour"]] <- getRedHairValues(
  phenotypeSources, phenotypeTables, correctionTable)

# Coronary artery disease
message("    Coronary artery disease...")
traitList[["Coronary artery disease"]] <- getCoronaryArteryDiseaseValues(
  phenotypeSources, phenotypeTables, correctionTable)

# Schizophrenia
message("    Schizophrenia...")
traitList[["Schizophrenia"]] <- getSchizophreniaValues(
  phenotypeSources, phenotypeTables, correctionTable)

message("    EduYears...")
traitList[["EduYears"]] <- getEduYears(
  phenotypeSources, phenotypeTables, correctionTable)

rm(phenotypeTables)
rm(correctionTable)
rm(phenotypeSources)
gc()

message("Completed processing individual traits!")
message("Plotting trait histograms...")
pdf(paste0(args$out, ".lifelines.debugHistograms.pdf"))
par(xpd = NA)
lapply(names(traitList), function(name) plotPhenotype(name, traitList[[name]]))
dev.off()

message(paste0("Debug histograms written to '", args$out, ".lifelines.debugHistograms.pdf", "'."))

# Combine to single table
lifelinesPhenotypeTable <- bind_rows(traitList, .id = "TRAIT") %>%
  ungroup() %>%
  select(PSEUDOIDEXT, AGE, SEX, VALUE, TRAIT)

ugliPhenotypeTable <- lifelinesPhenotypeTable %>%
  inner_join(gsaLinkageTable, by="PSEUDOIDEXT") %>%
  select(UGLI_ID, AGE, SEX, VALUE, TRAIT)

message("Writing output tables...")
write.table(lifelinesPhenotypeTable %>% rename(ID = PSEUDOIDEXT), 
            paste0(args$out, ".lifelines.dat"), row.names=F, col.names=T, quote=F, sep="\t")
message(paste0("Output written to '", args$out, ".lifelines.dat", "'."))

write.table(ugliPhenotypeTable %>% rename(ID = UGLI_ID), 
            paste0(args$out, ".ugli.dat"), row.names=F, col.names=T, quote=F, sep="\t")
message(paste0("Output written to '", args$out, ".ugli.dat", "'."))

message(paste0("DONE!"))