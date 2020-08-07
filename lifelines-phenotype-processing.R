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
default.phenotype_sources_map <- file.path(script.basename, "phenotype-sources-map.txt")
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
                    help='path to a file where the derived phenotypes are to be written to.')

##############################
# Define functions
##############################

loadPhenotypeTables <- function(filePath) {
  
  # phenotypeDataFrame <- fread(
  #   text = gsub("(?<!\r)(\n)", "\\n", read_file(filePath), perl = T), 
  #   header=T, quote="", sep="\t")

  phenotypeDataFrame <- fread(filePath,
    header=T, quote="", sep="\t", fill=T)
  
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
  
  return(phenotypeDataFrame)
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
  
  vmidAgeColumnMap <- sapply(vmidAgeNameMap, function(name) phenotypeSources$ColumnIdentifier[phenotypeSources$Name == name], USE.NAMES=T)
  
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
    inner_join(correctionTable, by = c("PSEUDOIDEXT", "VMID"))
  
  # Get follow up Mi
  columnIdentifierMiFollowUp <- phenotypeSources$ColumnIdentifier[phenotypeSources$Name == "Since last questionnaire: heart attack"]
  
  followUpMi <- phenotypeTables[[phenotypeSources$filePath[phenotypeSources$Name == "Since last questionnaire: heart attack"]]] %>%
    mutate(
      VALUE = case_when(get(columnIdentifierMiFollowUp) == "Yes" ~ 1L,
                        get(columnIdentifierMiFollowUp) == "No" ~ 0L,
                        TRUE ~ NA_integer_)) %>%
    filter(!is.na(VALUE)) %>%
    inner_join(correctionTable, by = c("PSEUDOIDEXT", "VMID"))
  
  # Get baseline intervention
  columnIdentifierInterventionBaseline <- phenotypeSources$ColumnIdentifier[phenotypeSources$Name == "Balloon angioplasty/bypass surgery"]
  
  baselineIntervention <- phenotypeTables[[phenotypeSources$filePath[phenotypeSources$Name == "Balloon angioplasty/bypass surgery"]]] %>%
    mutate(
      VALUE = case_when(get(columnIdentifierInterventionBaseline) == "Yes" ~ 1L,
                        get(columnIdentifierInterventionBaseline) == "No" ~ 0L,
                        TRUE ~ NA_integer_)) %>%
    filter(!is.na(VALUE)) %>%
    inner_join(correctionTable, by = c("PSEUDOIDEXT", "VMID"))
  
  columnIdentifierInterventionFollowUp <- phenotypeSources$ColumnIdentifier[phenotypeSources$Name == "Balloon angioplasty since last questionnaire"]
  
  FollowUpIntervention <- phenotypeTables[[phenotypeSources$filePath[phenotypeSources$Name == "Balloon angioplasty since last questionnaire"]]] %>%
    mutate(
      VALUE = case_when(get(columnIdentifierInterventionFollowUp) == "Yes" ~ 1L,
                        TRUE ~ NA_integer_)) %>%
    filter(!is.na(VALUE)) %>%
    inner_join(correctionTable, by = c("PSEUDOIDEXT", "VMID"))
  
  intermediateCad <- bind_rows(list(baselineMi, followUpMi, baselineIntervention, FollowUpIntervention))
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
                                         levels = VMID),
                  VALUE = get(columnIdentifier)) %>%
           filter(!is.na(VALUE)) %>%
           inner_join(correctionTable, by = c("PSEUDOIDEXT", "VMID")) %>%
           group_by(PSEUDOIDEXT) %>%
           slice(1) %>%
           select(PSEUDOIDEXT, AGE, SEX, VALUE))
}

getLatestValueFromRawPhenotypeTable <- function(phenotypeSources, phenotypeTables, correctionTable, name) {
  columnIdentifier <- phenotypeSources$ColumnIdentifier[phenotypeSources$Name == name]
  
  return(phenotypeTables[[phenotypeSources$filePath[phenotypeSources$Name == name]]] %>%
           mutate(VALUE = get(columnIdentifier),
                  VMID = if("VMID" %in% colnames(.)) VMID else getVmidFromEncountercode(ENCOUNTERCODE)) %>%
           filter(!is.na(VALUE)) %>%
           inner_join(correctionTable, by = c("PSEUDOIDEXT", "VMID")) %>%
           
           # Per PSEUDOIDEXT, get the row with the 'latest' VMID that is not NA
           group_by(PSEUDOIDEXT) %>%
           slice_max(as.numeric(VMID)) %>%
           select(PSEUDOIDEXT, AGE, SEX, VALUE))
}

getEstimatedGfr <- function(phenotypeSources, phenotypeTables, correctionTable) {
  creatinine <- getLatestValueFromRawPhenotypeTable(
    phenotypeSources, phenotypeTables, correctionTable, "Creatinine")

  estimatedGfr <- apply(creatinine, 1, function(row) {
    if (is.na(row["SEX"]) | (row["SEX"] != "Female" & row["SEX"] != "Male")) {
      return(NA)
    }
    
    return(estimateGfr(
      serum_creatine = as.numeric(row["VALUE"]), 
      age = as.numeric(row["age_bl1"]), 
      female = row["SEX"] == "Female", 
      black = FALSE))
  })
  
  return(creatinine %>%
    mutate(VALUE = estimatedGfr) %>%
    filter(is.na(VALUE)))
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

##############################
# Generate phenotypes
##############################

args <- parser$parse_args(commandArgs(trailingOnly = TRUE))

args <- parser$parse_args(c(
  "--gsa-linkage-file", "/groups/umcg-lifelines/tmp01/releases/gsa_linkage_files/v1/gsa_linkage_file.dat",
  "--out", "/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/data/lifelines/processed/UGLI.pgs.phenotypes.dat",
  "--phenotype-source-map", "/home/umcg-rwarmerdam/pgs_based_mixup_correction/scripts/r-scripts/pgs_based_sample_mix-up_correction/phenotype-source-map.txt"))

# Load GSA linkage file
gsaLinkageFile <- fread(args$gsa_linkage_file, sep="\t", header=T)

# Collect required columns
phenotypeSources <- fread(args$phenotype_source_map, sep="\t", header=T) %>%
  mutate(filePath = file.path(Path, FileName))

# For every unique file, load it, and get the respective columns.
phenotypeTables <- sapply(unique(phenotypeSources$filePath),
       loadPhenotypeTables, USE.NAMES = TRUE, simplify = F)

correctionTable <- getCorrectionTable(phenotypeSources, phenotypeTables)

traitList <- list()

# Height
traitList[["Height"]] <- getLatestValueFromRawPhenotypeTable(
  phenotypeSources, phenotypeTables, correctionTable, "Height")

# BMI
traitList[["BMI"]] <- getLatestValueFromRawPhenotypeTable(
  phenotypeSources, phenotypeTables, correctionTable, "BMI")

# Diastolic blood pressure
traitList[["Diastolic blood pressure"]] <- getBloodPressureValues(
  phenotypeSources, phenotypeTables, correctionTable, "Diastolic blood pressure")

# Systolic blood pressure
traitList[["Systolic blood pressure"]] <- getBloodPressureValues(
  phenotypeSources, phenotypeTables, correctionTable, "Systolic blood pressure")

# Estimated GFR
traitList[["eGFR"]] <- getEstimatedGfr(
  phenotypeSources, phenotypeTables, correctionTable)

# Basophil concentration
traitList[["Basophilic Granulocytes"]] <- getLatestValueFromRawPhenotypeTable(
  phenotypeSources, phenotypeTables, correctionTable, "Basophilic Granulocytes")

# Eosinophil concentration
traitList[["Eosinophil concentration"]] <- getLatestValueFromRawPhenotypeTable(
  phenotypeSources, phenotypeTables, correctionTable, "Eosinophil concentration")

# HbA1c concentration
traitList[["HbA1c"]] <- getLatestValueFromRawPhenotypeTable(
  phenotypeSources, phenotypeTables, correctionTable, "HbA1c")

# LDL-cholesterol
traitList[["LDL cholesterol"]] <- getLatestValueFromRawPhenotypeTable(
  phenotypeSources, phenotypeTables, correctionTable, "LDL Cholesterol")

# Log-transformed triglyceride concentration
traitList[["triglyceride concentration (log-transformed)"]] <- getLogTransformedTriglycerideConcentration(
  phenotypeSources, phenotypeTables, correctionTable)

# Square root of HDL-cholesterol
traitList[["HDL cholesterol (square root)"]] <- getSquareRootOfHdlCholesterol(
  phenotypeSources, phenotypeTables, correctionTable)

# Total cholesterol concentration
traitList[["Total cholesterol"]] <- getLatestValueFromRawPhenotypeTable(
  phenotypeSources, phenotypeTables, correctionTable, "Cholesterol")

# Hematocrit concentration
traitList[["Hematocrit concentration"]] <- getLatestValueFromRawPhenotypeTable(
  phenotypeSources, phenotypeTables, correctionTable, "Hematocrit")

# Hemoglobin concentration
traitList[["Hemoglobin concentration"]] <- getLatestValueFromRawPhenotypeTable(
  phenotypeSources, phenotypeTables, correctionTable, "Hemoglobin")

# Lymphocyte concentration
traitList[["Lymphocyte concentration"]] <- getLatestValueFromRawPhenotypeTable(
  phenotypeSources, phenotypeTables, correctionTable, "Lymphocytes")

# Monocyte concentration
traitList[["Monocyte concentration"]] <- getLatestValueFromRawPhenotypeTable(
  phenotypeSources, phenotypeTables, correctionTable, "Monocytes")

# Erythrocyte concentration
traitList[["Erythrocyte concentration"]] <- getLatestValueFromRawPhenotypeTable(
  phenotypeSources, phenotypeTables, correctionTable, "Erythrocytes")

# Neutrophil concentration
traitList[["Neutrophil concentration"]] <- getLatestValueFromRawPhenotypeTable(
  phenotypeSources, phenotypeTables, correctionTable, "Neutrophil Granulocytes")

# Thrombocyte concentration
traitList[["Thrombocyte concentration"]] <- getLatestValueFromRawPhenotypeTable(
  phenotypeSources, phenotypeTables, correctionTable, "Thrombocytes")

# Type-2 diabetes
traitList[["Type 2 diabetes"]] <- getDerivedBinaryValues(
  phenotypeSources, phenotypeTables, correctionTable, "Type 2 diabetes", "Second assessment (2A)")

# Type-1 diabetes
traitList[["Type 1 diabetes"]] <- getDerivedBinaryValues(
  phenotypeSources, phenotypeTables, correctionTable, "Type 1 diabetes", "Second assessment (2A)")

# Asthma
traitList[["Asthma"]] <- getDerivedBinaryValues(
  phenotypeSources, phenotypeTables, correctionTable, "Asthma", "Baseline assessment (1A)")

# Hair colour (level of blonde (black - blonde))
traitList[["Blondeness of hair"]] <- getBlondenessOfHair(
  phenotypeSources, phenotypeTables, correctionTable)

# Hair colour (red vs. a level of blonde (black - blonde))
traitList[["Red hair colour"]] <- getRedHairValues(
  phenotypeSources, phenotypeTables, correctionTable)

# Coronary artery disease
traitList[["Coronary artery disease"]] <- getCoronaryArteryDiseaseValues(
  phenotypeSources, phenotypeTables, correctionTable)

# Schizophrenia
traitList[["Schizophrenia"]] <- getSchizophreniaValues(
  phenotypeSources, phenotypeTables, correctionTable)

bind_rows(traitList, .id = "TRAIT")