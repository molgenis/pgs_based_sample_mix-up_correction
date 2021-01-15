#!/usr/bin/env Rscript

# Library
library(argparse)

# Data transformation
label.log_egfr <- "logegfr"
label.log_triglycerides <- "logtgl"
label.sqrt_hdl <- "sqrthdc"

col.ifabsent <- function(phenotypes, col.name, values) {
  if (col.name %in% colnames(phenotypes)) {
    print(paste0("'", col.name, "' already present, skipping..."))
  } else {
    phenotypes[,col.name] <- values
    print(paste0("Added '", col.name, "'"))
  }
}

parser <- ArgumentParser(description='add derivatives to the phenotypes')
parser$add_argument('--path',
                    help='path to phenotype table')
parser$add_argument('--out',
                    help='output path')

args <- parser$parse_args(commandArgs(trailingOnly = TRUE))

phenotypes <- read.table(args$path, header=T, sep="\t", row.names=1, quote="")

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

estimatedGFR <- apply(phenotypes, 1, function(row) {
  if (is.na(row["geslacht"]) | (row["geslacht"] != "Female" & row["geslacht"] != "Male")) {
    return(NA)
  }
  
  return(estimateGFR(as.numeric(row["bkr"]), as.numeric(row["age_bl1"]), row["geslacht"] == "Female", FALSE))
  })

col.ifabsent(phenotypes, label.log_egfr, log(estimatedGFR))

col.ifabsent(phenotypes, label.log_triglycerides, log(phenotypes$tgl))

col.ifabsent(phenotypes, label.sqrt_hdl, sqrt(phenotypes$hdc))

# Add diseases



write.table(phenotypes, args$out, col.names=T, sep="\t", row.names=T, quote=F)

