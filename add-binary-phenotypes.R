#!/usr/bin/env Rscript

# Library
library(argparse)

# Labels
#label.schizophrenia <- "health72i9.recode"
#label.bipolar <- "health72i8.recode"
label.t2d <- "health17.recode.t2"
label.t1d <- "health17.recode.t1"
#label.depression <- "health72i3.recode"
#label.rheumatiod <- "health_72g3a.recode"
label.asthma <- "health12.recode"
#label.celiac <- "health72c7.recode"
#label.alzheimers <- "health_72o.recode"
# label.sle <- "csq21.recode"
label.pollen <- "health15d.recode"

col.ifabsent <- function(phenotypes, col.name, values) {
  if (col.name %in% colnames(phenotypes)) {
    print(paste0("'", col.name, "' already present, skipping..."))
  } else {
    phenotypes[,col.name] <- values
    print(paste0("Added '", col.name, "'"))
  }
  return(phenotypes)
}

parser <- ArgumentParser(description='add derivatives to the phenotypes')
parser$add_argument('--path',
                    help='path to phenotype table')
parser$add_argument('--out',
                    help='output path')

args <- parser$parse_args(commandArgs(trailingOnly = TRUE))

#phenotypes <- read.table("LLDeep_GoNL_samples_V01.txt", header=T, sep="\t", quote="")
phenotypes <- read.table(args$path, header=T, sep="\t", quote="")

# ASTHMA
# ======
# Get true logical value of people that have asthma and are
# diagnosed by a doctor.
asthma.bool <- phenotypes["health12a"] == "Yes"
# People that don't have asthma should get false instead
# of NA. NAs in the asthma field should be maintained
asthma.bool[is.na(asthma.bool) & !is.na(phenotypes["health12"])] <- FALSE
# Map false to 0 and true to 1.
asthma.01 <- as.integer(asthma.bool)
phenotypes <- col.ifabsent(phenotypes, label.asthma, asthma.01)

# DIABETES TYPE 2
# ===============
# Get true where an individual has diabetes and an
# individual has type 2 diabetes. 
# NAs where diabetes (health17a) is NA
t2d <- as.integer(as.logical(
  phenotypes["health17a"] == "Yes" 
  & phenotypes["health17b"] == "Type 2 (adult-onset diabetes, later in life)"))
phenotypes <- col.ifabsent(phenotypes, label.t2d, t2d)

# DIABETES TYPE 1
# ===============
# Get true where an individual has diabetes and an
# individual has type 1 diabetes. 
# NAs where diabetes (health17a) is NA
t1d <- as.integer(as.logical(
  phenotypes["health17a"] == "Yes" 
  & phenotypes["health17b"] == "Type 1 (juvenile diabetes, from childhood)"))
phenotypes <- col.ifabsent(phenotypes, label.t1d, t1d)

# SLE / LUPUS
# ===========

# as.integer(as.logical(
#   phenotypes["csq21"] == "Yes" 
#   & phenotypes["csq21a"] == "Positive (lupus/SLE diagnosed)"))

# POLLEN
# ======

pollen <- as.integer(
  !is.na(phenotypes["health15d"]) 
  & phenotypes["health15d"] == "Yes")
phenotypes <- col.ifabsent(phenotypes, label.pollen, pollen)

write.table(phenotypes, args$out, row.names=F, col.names=T, quote=F, sep="\t")