#!/usr/bin/env Rscript

library(argparse)

parser <- ArgumentParser(description='')
parser$add_argument('--phenotypes_path',
                    help='phenotype table')
parser$add_argument('--gwas_map_path',
                    help='file with phenotype ids')
parser$add_argument('--to_correct_for', nargs="+",
                    help='traits / data to correct phenotypes for')
parser$add_argument('--out', help='output file')

command.args <- commandArgs(trailingOnly = TRUE)

if (length(command.args) > 0) {
  # Parse arguments
  args <- parser$parse_args(command.args)
  
  phenotypes_path <- args$phenotypes_path
  
  phenotypes <- read.table(phenotypes_path, header=T, quote="", sep="\t")

  gwas_map_path <- args$gwas_map_path

  gwas_map <- read.table(gwas_map_path, header=F, quote="", sep=",", col.names = c("id", "label"), stringsAsFactors = F, row.names = "id")
  to_correct_for <- args$to_correct_for

  out <- args$out
  
} else {
  height <- read.table()
  sex <- read.table()
  age <- read.table()
  
  phenotypes <- merge(height, sex, by="row.names")
  rownames(phenotypes) <- phenotypes$Row.names
  phenotypes$Row.names <- NULL
  
  phenotypes <- merge(phenotypes, age, by="row.names")
  rownames(phenotypes) <- phenotypes$Row.names
  phenotypes$Row.names <- NULL
  
  gwas_map <- data.frame(id = c("health17.recode.t2"), 
                         label = c("t2d"), 
                         stringsAsFactors = F)
  
  to_correct_for <- c("geslacht", "age_bl1")
  
  out <- "~/Documents/out"
}

rownames(gwas_map) <- gwas_map$id

correct <- function(df) {
  colnames(df) <- c("y", colnames(df)[2:ncol(df)])
  model.2 <- lm(y ~ . + .^2, df)
  print(summary(model.2))
  return(resid(model.2))
}

for (id in rownames(gwas_map)) {
  id.adj <- paste0(id, ".adj")
  print(id)
  phenotypes.adj <- correct(phenotypes[, c(id, to_correct_for)])
  names.adj <- names(phenotypes.adj)
  phenotypes[id.adj] <- NA
  phenotypes[names.adj, id.adj] <- phenotypes.adj
  gwas_map[id, "id"] <- id.adj
  print(paste0("adjusted ",id))
}

write.table(gwas_map[c("label", "id")], paste0(out, "gwasToTrait.csv"), row.names=F, sep=",", col.names=F, quote=F)
write.table(phenotypes, paste0(out, "phenotypes.tsv"), row.names=F, sep="\t", col.names=T, quote=F)

