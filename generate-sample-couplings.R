#!/usr/bin/env Rscript

#############################################################
## c.a.warmerdam@umcg.nl - September 2020
## Generates sample coupling file.
#############################################################

##############################
# Load libraries
##############################
library(tidyverse)
library(argparse)
library(data.table)

##############################
# Define argument parser
##############################
parser <- ArgumentParser(description='')

parser$add_argument('--mix-up-percentage', required=FALSE, type="double",
                    help=paste0('introduce mix-ups in link file',
                                'and phenotype sample ids in the second column'))
parser$add_argument('--sample-count', required=FALSE, type="integer", default="5120",
                    help=paste0('number of samples to include in the coupling file'))
parser$add_argument('--out', help='path to output prefix', required=T)

inputGroup <- parser$add_mutually_exclusive_group(required=T)

inputGroup$add_argument('--sample-coupling-file',
                        help=paste0('file containing genotype sample ids in the first column',
                                    'and phenotype sample ids in the second column'))
inputGroup$add_argument('--phenotypes-file',
                        help='path to a tab-delimited file holding all processed phenotype data.')

##############################
# Define functions
##############################
permute <- function(linkTable, nMixUpsToIntroduce) {
  linkTable$permuted <- linkTable$geno
  
  stopifnot("Cannot mix-up a single sample!" = 1 < nMixUpsToIntroduce)

  # Shuffle all samples so the ones in which mix-ups are introduced are random
  shuffledSamples <- sample(linkTable$permuted)
  
  # Create a vector of booleans indicating whether or not a sample has been mixed up.
  swappedSamples <- rep(F, length(shuffledSamples))
  
  # Perform the code below for every mix-up to introduce, from n to 1 (inclusive).
  for (i in rev(2:(nMixUpsToIntroduce))) {
    # Should the current sample (i) be mixed up?
    # Swap the sample if either: the sample is not yet swapped, or
    # the i sample is the second-last sample to swap and the last sample has to be swapped still.
    mustSwapCurrent <- (!swappedSamples[i] | (i == 2 & !swappedSamples[1]))

    j <- i
    
    # Set the sample to swap the current one with.
    if (mustSwapCurrent) {
      # If the current one must be swapped, the sample to swap with
      # should be the one of the '1'st to the 'i'th (exclusive).
      j <- sample.int(i - 1, size = 1)
    } else {
      # If the current one does not have to be swapped, 
      # the sample to swap with should be the one of the '1'st to the 'i'th (inclusive).
      j <- sample.int(i, size = 1);
    }
    
    # If j is equal to i, the 'i'th sample is already swapped.
    if (j != i) {
      
      # Get the samples to swap
      a <- shuffledSamples[i]
      b <- shuffledSamples[j]

      # Replace the samples in the map.
      aLoc <- linkTable$permuted == a
      bLoc <- linkTable$permuted == b
      linkTable$permuted[aLoc] <- b
      linkTable$permuted[bLoc] <- a
      
      # Swap the samples in the shuffled sample list.
      shuffledSamples[c(i, j)] <- shuffledSamples[c(j, i)]
      
      # Set the samples to shuffled.
      swappedSamples[i] <- T
      swappedSamples[j] <- T
    }
  }
  
  return(linkTable)
}

##############################
# Run
##############################
args <- parser$parse_args(commandArgs(trailingOnly = TRUE))

# Get the file with 
sampleCouplingFilePath <- args$sample_coupling_file

# Get the 
phenotypesFilePath <- args$phenotypes_file

# Get the percentage of mix-ups to introduce / induce
mixUpPercentage <- args$mix_up_percentage

# Number of samples to output
numberOfSamples <- args$sample_count

out <- args$out

link <- NULL

if (!is.null(sampleCouplingFilePath)) {
  link <- fread(sampleCouplingFilePath, stringsAsFactors=F,
                col.names = c("geno", "pheno")) %>%
    filter(pheno %in% unique(phenotypesTable$ID))
}

if (!is.null(phenotypesFilePath)) {
  phenotypesTable <- fread(phenotypesFilePath, header=T, quote="", sep="\t",
                           col.names = c("ID", "AGE", "SEX", "VALUE", "TRAIT")) %>%
    mutate(SEX = factor(SEX, levels = c("Female", "Male"))) %>%
    group_by(ID) %>%
    filter(!any(AGE < 18)) %>%
    ungroup()
  
    link <- data.frame(geno = unique(phenotypesTable$ID), pheno = unique(phenotypesTable$ID))
}

if (numberOfSamples & numberOfSamples < nrow(link)) {
  link <- link %>%
    slice_sample(n = numberOfSamples)
}

if (mixUpPercentage & mixUpPercentage > 0 & mixUpPercentage <= 10) {
  nMixUpsToIntroduce <- round(
    nrow(link) / 100 * mixUpPercentage);
  
  nMixUpsToIntroduce <- max(nMixUpsToIntroduce, 2)
  
  permutedLink <- permute(link, nMixUpsToIntroduce)
  
  numberOfPermutedSamples <- sum(permutedLink$geno != permutedLink$permuted)
  
  message(paste0("Mix-ups introduced: ", numberOfPermutedSamples, "/", nrow(permutedLink), 
                 " (", numberOfPermutedSamples / nrow(permutedLink) * 100, "%)"))
  
  permutedLink <- rename(permutedLink, original = geno, geno = permuted) %>%
    select(pheno, geno, original)
  
  write.table(permutedLink %>% filter(original != geno), 
              paste0(out, ".perm_", numberOfPermutedSamples, "mixUps.txt"),
              row.names = F, col.names = T, quote = F, sep = "\t")
  
  write.table(permutedLink, 
              paste0(out, ".perm_", numberOfSamples, "samples_", numberOfPermutedSamples, "mixUps.txt"), 
              row.names = F, col.names = T, quote = F, sep = "\t")
}

write.table(link, paste0(out, ".", numberOfSamples, "samples.txt"), 
            row.names = F, col.names = T, quote = F, sep = "\t")
