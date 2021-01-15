#!/usr/bin/env Rscript

#############################################################
## c.a.warmerdam@umcg.nl - November 2020
## Assess difference in AUC between NB methods per trait
#############################################################

##############################
# Load libraries
##############################
library(tidyverse)
library(data.table)
library(ggpubr)
library(Hmisc)
library(ggforce)

##############################
# Define argument parser
##############################

parser <- ArgumentParser(description='')
parser$add_argument('files', nargs = "+",
                    help='Output statistic files (\'outputStatisticsPerTrait.tsv\') to include in likelihood model comparison.')

##############################
# Run
##############################

old <- theme_set(theme_classic())
theme_update(line = element_line(
  colour = "black", size = (0.5 * ggplot2::.pt * 72.27/96), 
  linetype = 1, lineend = "butt", arrow = F, inherit.blank = T),
  strip.background = element_rect(colour = NA, fill = NA),
  axis.line = element_line(
    colour = "#595A5C", size = (0.5 * ggplot2::.pt * 72.27/96), 
    linetype = 1, lineend = "butt", arrow = F, inherit.blank = T),
  axis.ticks = element_line(
    colour = "#595A5C", size = (0.5 * ggplot2::.pt * 72.27/96), 
    linetype = 1, lineend = "butt", arrow = F, inherit.blank = T)
)

scientific_10 <- function(x) {
  gsub("e", " %*% 10^", scales::scientific_format()(x))
}

args <- parser$parse_args(commandArgs(trailingOnly = TRUE))
runs <- tibble(fullPath = args$files, label = seq_along(args$files))

# Define function for reading in output statistics per trait
lapply_read_csv_bind_rows <- function(files, labels) {
  lapply(seq_along(files), function(filePathIndex) {
    read.table(files[filePathIndex], sep = "\t") %>%
      mutate(run = labels[filePathIndex])
    }) %>% bind_rows()
}

# Load all output statistics from parameter sweeps
resultsTable <- lapply_read_csv_bind_rows(runs$fullPath, runs$label) %>%
  filter(included == "Y") %>%
  select(trait, traitDataType, confinedAuc, matrixWideAuc, confinedAucOnScaledLlr, pValue, naiveBayesMethod, samplesPerNaiveBayesBin, run) %>%
  mutate(label = paste0(naiveBayesMethod, ".", samplesPerNaiveBayesBin))

# Create summary of results with means, sd, min and max
resultsTableSummary <- resultsTable %>%
  group_by(trait, label) %>%
  summarise(confinedAucMean = mean(confinedAuc), confinedAucSd = sd(confinedAuc), 
            confinedScaledAucMean = mean(confinedAucOnScaledLlr), confinedScaledAucSd = sd(confinedAucOnScaledLlr), 
            confinedScaledAucMin = min(confinedAucOnScaledLlr), confinedScaledAucMax = max(confinedAucOnScaledLlr),
            matrixWideAucMean = mean(matrixWideAuc), matrixWideAucSd = sd(matrixWideAuc),
            pValue = mean(pValue), pValue = sd(pValue),
            n = n())

# resultsTableTesting <- resultsTableSummary %>%
#   full_join(resultsTableSummary, by = "trait") %>%
#   group_by(trait, label.x, label.y) %>%
#   mutate(pVal = t.test(
#     resultsTable[resultsTable$trait == trait & resultsTable$label == label.x, "confinedAucOnScaledLlr"], 
#     resultsTable[resultsTable$trait == trait & resultsTable$label == label.y, "confinedAucOnScaledLlr"])$p.value) %>%
#   select(trait, label.x, label.y, pVal)

getAnovaResults <- function(aov.res) {
  multiplePairwise <- TukeyHSD(aov.res)$label

  gaussianPairwiseNames <- c("gaussian.NA-efi-discretization.20", "gaussian.NA-efi-discretization.30",
                             "gaussian.NA-efi-discretization.50", "gaussian.NA-efi-discretization.80",
                             "gaussian.NA-ewi-discretization.20", "gaussian.NA-ewi-discretization.30",
                             "gaussian.NA-ewi-discretization.50", "gaussian.NA-ewi-discretization.80")

  # For the parameters with the highest mean AUC, we want to assess significant difference compared to other traits.
  gaussianPairwise <- as_tibble(multiplePairwise[gaussianPairwiseNames,]) %>%
    mutate(label = c("efi-discretization.20", "efi-discretization.30",
                     "efi-discretization.50", "efi-discretization.80",
                     "ewi-discretization.20", "ewi-discretization.30",
                     "ewi-discretization.50", "ewi-discretization.80")) %>%
    select(c(label, "diff" = "diff", "pAdjusted" = "p adj"))
}

anovaTestResultsPerTrait <- resultsTable %>%
  group_by(trait) %>%
  summarise(aovPval = summary(aov(confinedAucOnScaledLlr ~ label))[[1]][["Pr(>F)"]][1])

postHocAovResults <- mapply(function(thisTrait) {
  return(getAnovaResults(aov(confinedAucOnScaledLlr ~ label, 
                      data = resultsTable %>% filter(trait == thisTrait))) %>%
    mutate(trait = thisTrait))
  }, anovaTestResultsPerTrait$trait, SIMPLIFY = F) %>% bind_rows() %>%
  mutate(significanceLabel = case_when(
    pAdjusted <= 0.0001 ~ "****",
    pAdjusted <= 0.001 ~ "***",
    pAdjusted <= 0.01 ~ "**",
    pAdjusted <= 0.05 ~ "*",
    pAdjusted > 0.05 ~ ""))

gaussianParameterColours = c(
  "#009445",
  "#3980A6", "#19709E", "#05486C", "#023652",
  "#FF774F", "#F94E1C", "#AB2600", "#811D00")

names(gaussianParameterColours) <- c(
  "gaussian.NA", 
  paste0("efi-discretization.", c(20, 30, 50, 80)), 
  paste0("ewi-discretization.", c(20, 30, 50, 80)))

ggplot(resultsTableSummary, aes(x = trait, y = confinedScaledAucMean, group = label, colour = label)) +
  geom_point(position=position_dodge(0.5)) +
  geom_errorbar(aes(ymin = confinedScaledAucMean - confinedScaledAucSd,
                    ymax = confinedScaledAucMean + confinedScaledAucSd),
                position=position_dodge(0.5)) +
  scale_colour_manual(values = gaussianParameterColours) +
  scale_fill_manual(values = gaussianParameterColours) +
  theme(axis.text.x = element_text(angle = 90))

pdf(paste0("output/figures/likelihoodClassifierComparison_", format(Sys.Date(), "%Y%m%d"), ".pdf"), 
    useDingbats = FALSE, width = 8.5, height = 11)

par(xpd = NA)

ggplot(resultsTable, 
       aes(x = label, y = confinedAucOnScaledLlr, colour = label)) +
  geom_violin() +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), 
               geom="pointrange") +
  scale_colour_manual(values = gaussianParameterColours) +
  scale_fill_manual(values = gaussianParameterColours) +
  # scale_x_discrete(limits = sub(pattern = ".", replacement = ", ", x = unique(resultsTable$label)))
  geom_text(
    data = resultsTableSummary %>% left_join(postHocAovResults, by = c("trait", "label")),
    aes(y = confinedScaledAucMax, label = significanceLabel), 
    vjust = -0.12, show.legend = NA) +
  geom_label(
    data = resultsTableSummary %>% group_by(trait) %>% summarise(minAuc = min(confinedScaledAucMin)) %>% inner_join(anovaTestResultsPerTrait, by = "trait"),
    parse = T, aes(label = paste0("Anova:~italic(p)==", scientific_10(aovPval)), y = minAuc, x = 9.5), 
    vjust = 0.86, hjust = "inward", show.legend = NA, inherit.aes = F, size = 0.5 * 72.27/96 * 8, 
    fill = "#E9E9E8",   label.padding = unit(0.1, "lines"), label.r = unit(0.0, "lines"),
    label.size = 0.0) +
  scale_y_continuous(expand = expansion(mult = c(0.10, 0.12), add = c(0.024, 0))) +
  theme(axis.text.x = element_text(angle = 90), legend.position = "null") +
  ylab("Area under ROC curve") +
  xlab("Likelihood model") +
  labs(title = "Comparison of likelihood models", subtitle = "A: 1-15 / 25 traits") +
  facet_wrap_paginate(~ trait, scales = "free_y", ncol = 3, nrow = 5, page = 1)

ggplot(resultsTable, 
       aes(x = label, y = confinedAucOnScaledLlr, colour = label)) +
  geom_violin() +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1),
               geom="pointrange") +
  scale_colour_manual(values = gaussianParameterColours) +
  scale_fill_manual(values = gaussianParameterColours) +
  geom_text(
    data = resultsTableSummary %>% left_join(postHocAovResults, by = c("trait", "label")),
    aes(y = confinedScaledAucMax, label = significanceLabel), 
    vjust = -0.12, show.legend = NA) +
  geom_label(
    data = resultsTableSummary %>% group_by(trait) %>% summarise(minAuc = min(confinedScaledAucMin)) %>% inner_join(anovaTestResultsPerTrait, by = "trait"),
    parse = T, aes(label = paste0("Anova:~italic(p)==", scientific_10(aovPval)), y = minAuc, x = 9.5), 
    vjust = 0.86, hjust = "inward", show.legend = NA, inherit.aes = F, size = 0.5 * 72.27/96 * 8, 
    fill = "#E9E9E8",   label.padding = unit(0.1, "lines"), label.r = unit(0.0, "lines"),
    label.size = 0.0) +
  scale_y_continuous(expand = expansion(mult = c(0.10, 0.12), add = c(0.024, 0))) +
  theme(axis.text.x = element_text(angle = 90), legend.position = "null") +
  ylab("Area under ROC curve") +
  xlab("Likelihood model") +
  labs(title = "Comparison of likelihood models", subtitle = "B: 16-25 / 25 traits") +
  facet_wrap_paginate(~ trait, scales = "free_y", ncol = 3, nrow = 5, page = 2)

dev.off()
