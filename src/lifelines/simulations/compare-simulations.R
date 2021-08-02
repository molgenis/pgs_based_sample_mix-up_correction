#!/usr/bin/env Rscript

#############################################################
## c.a.warmerdam@umcg.nl - July 2021
## Assess difference in AUC with various numbers of phenotypes
#############################################################

##############################
# Load libraries
##############################
library(tidyverse)
library(data.table)
library(ggforce)
library(Hmisc)

##############################
# Constants
##############################


datasets <- c(
  "simulationDataset_nTraits10_evf0.5_20210719", "simulationDataset_nTraits25_evf0.5_20210719",
  "simulationDataset_nTraits50_evf0.5_20210719", "simulationDataset_nTraits75_evf0.5_20210719",
  "simulationDataset_nTraits100_evf0.5_20210719", "simulationDataset_nTraits10_evf1_20210703", 
  "simulationDataset_nTraits25_evf1.5_20210715", "simulationDataset_nTraits50_evf2_20210703",
  "simulationDataset_nTraits100_evf1_20210703", "simulationDataset_nTraits10_evf1.5_20210715",
  "simulationDataset_nTraits25_evf2_20210703", "simulationDataset_nTraits75_evf1_20210703",
  "simulationDataset_nTraits100_evf1.5_20210715", "simulationDataset_nTraits10_evf2_20210703", 
  "simulationDataset_nTraits50_evf1_20210703", "simulationDataset_nTraits75_evf1.5_20210715",
  "simulationDataset_nTraits100_evf2_20210703", "simulationDataset_nTraits25_evf1_20210703", 
  "simulationDataset_nTraits50_evf1.5_20210715", "simulationDataset_nTraits75_evf2_20210703")

iterations <- c("20210709", "20210709-2", "20210709-3", "20210709-4", "20210709-5")

resultsNoTraits <- tibble(confinedAucLlrSexCheck = 0.7415957, nTraits = 0, evf = c(0.5, 1, 1.5, 2))

dataPath <- "~/Documents/projects/pgs_based_mixup_correction-ugli/data/simulated/"

outputPath <- "~/Documents/projects/pgs_based_mixup_correction-ugli/output/simulations"

## Plotting defaults
old <- theme_set(theme_classic())
theme_update(line = element_line(
  colour = "black", size = (0.5 / (ggplot2::.pt * 72.27/96)),
  linetype = 1, lineend = "butt", arrow = F, inherit.blank = T),
  strip.background = element_rect(colour = NA, fill = NA),
  axis.line = element_line(
    colour = "#595A5C", size = (0.5 / (ggplot2::.pt * 72.27/96)),
    linetype = 1, lineend = "butt", arrow = F, inherit.blank = T),
  axis.ticks = element_line(
    colour = "#595A5C", size = (0.5 / (ggplot2::.pt * 72.27/96)),
    linetype = 1, lineend = "butt", arrow = F, inherit.blank = T)
)

##############################
# Run
##############################

mean_whiskers <- function(x) {
  return(c("y" = mean(x), "ymin" = min(x), "ymax" = max(x)))
}

lapply_read_csv_bind_rows <- function(files, labels) {
  lapply(seq_along(files), function(filePathIndex) {
    read.table(files[filePathIndex], sep = "\t") %>%
      mutate(run = labels[filePathIndex])
  }) %>% bind_rows()
}

fileTibble <- expand_grid(datasets, iterations) %>%
  mutate(files = file.path(
    outputPath,
    datasets, "split_prediction", iterations, "2.predictAuc/overallOutputStatistics.tsv"))
files <- fileTibble$files

runs <- tibble(fullPath = files, label = seq_along(files), iterations = fileTibble$iterations) %>%
  filter(file.exists(files))

matches <- str_match(runs$fullPath, "simulationDataset_nTraits(10|25|50|75|100)_evf([:digit:](?:\\.[:digit:]+)?)_[:digit:]{8}")
runs$dataSetLab <- matches[,1]
runs$nTraits <- as.numeric(matches[,2])
runs$evf <- as.numeric(matches[,3])

runsCorr <- runs %>% 
  mutate(actualCorrelationsPath = file.path(
    "/Users/cawarmerdam/Documents/projects/pgs_based_mixup_correction-ugli/data/simulated", 
    runs$dataSetLab, "actualCorrelations.txt")) %>%
  select(actualCorrelationsPath, dataSetLab) %>% distinct()

resultsTable <- lapply_read_csv_bind_rows(runs$fullPath, runs$label) %>%
  mutate(confinedAuc = case_when(confinedAuc >= 0.5 ~ confinedAuc, TRUE ~ 1 - confinedAuc),
         confinedAucLlrSexCheck = case_when(confinedAucLlrSexCheck >= 0.5 ~ confinedAucLlrSexCheck, TRUE ~ 1 - confinedAucLlrSexCheck)) %>%
  inner_join(runs, by = c("run" = "label")) %>%
  mutate(run = as.character(run)) %>%
  bind_rows(resultsNoTraits)

simulationNameTable <- c(
  "SIM1" = "1st simulated version",
  "SIM2" = "2nd simulated version",
  "SIM3" = "3rd simulated version",
  "SIM4" = "4th simulated version")

corrTable <- lapply(seq_along(runsCorr$actualCorrelationsPath), function(filePathIndex) {
    actualCorrelations <- as.data.frame.table(
      as.matrix(read.table(runsCorr$actualCorrelationsPath[filePathIndex], 
                           sep = "\t", check.names = F)), 
      responseName = "correlations") %>%
      separate(Var1, c("rowType", "rowName", "rowSimName"), sep = "_") %>%
      separate(Var2, c("colType", "colName", "colSimName"), sep = "_") %>%
      filter(rowName == colName, rowSimName == colSimName, rowType == "PGS", colType == "VALUE") %>%
      mutate(fullName = paste0(rowName, " (", simulationNameTable[rowSimName], ")")) %>%
      dplyr::select(rowSimName, fullName, correlations)
    colnames(actualCorrelations) <- c("rowSimName", "fullName", as.character(runsCorr$dataSetLab[filePathIndex]))
    return(actualCorrelations)
  }) %>% reduce(full_join, by = c("fullName", "rowSimName")) %>%
  arrange(rowSimName)

resultsTableWide <- resultsTable %>% 
  dplyr::select(nTraits, evf, dataSetLab, iterations, confinedAuc, confinedAucLlrSexCheck) %>%
  filter(!is.na(iterations)) %>%
  pivot_wider(id_cols = c(nTraits, evf, dataSetLab), 
              names_from = iterations, 
              names_glue = "{iterations}_{.value}", 
              values_from = c(confinedAuc, confinedAucLlrSexCheck)) %>%
  rowwise() %>% 
  mutate(`Mean confinedAuc` = mean(c_across(cols = ends_with("confinedAuc"))),
         `Mean confinedAucLlrSexCheck` = mean(c_across(cols = ends_with("confinedAucLlrSexCheck"))))

corrTableWide <- corrTable %>% 
  pivot_longer(c(-rowSimName, -fullName), names_to = "dataSetLab", values_to = "value") %>%
  inner_join(resultsTableWide, by = "dataSetLab") %>%
  pivot_wider(id_cols = -rowSimName, names_from = "fullName", values_from = "value") %>%
  dplyr::select(nTraits, evf, ends_with("confinedAuc"), ends_with("confinedAucLlrSexCheck"), corrTable$fullName) %>%
  arrange(nTraits, evf)

write.table(corrTableWide, 
            paste0("output/tables/simulationsComparison_", format(Sys.Date(), "%Y%m%d"), ".txt"), 
    quote = F, row.names = F, col.names = T, sep = "\t")

pdf(paste0("output/figures/simulationsComparison_", format(Sys.Date(), "%Y%m%d"), ".pdf"), 
    useDingbats = FALSE, width = 6, height = 4.2)

par(xpd = NA)

colorFun <- colorRampPalette(c("#D19D00", "#0F358C"))

ggplot(resultsTable, aes(x = nTraits, y = confinedAucLlrSexCheck, group = factor(evf), colour = factor(evf), fill = factor(evf))) +
  geom_hline(yintercept = 1, colour = "#595A5C", linetype = "dashed", size = (0.5 / (ggplot2::.pt * 72.27/96))) +
  stat_summary(fun=mean, 
               geom="line") +
  stat_summary(fun.data=mean_whiskers,
               geom="errorbar", width=1.2) +
  stat_summary(fun=mean, geom="point") +
  #geom_point(fill = "#ffffff", shape = 2, alpha = 0.5) +
  scale_colour_manual(name = "Explained variance factor", 
                      values = colorFun(4)) +
  geom_point(data = tibble(confinedAucLlrSexCheck = 0.90, nTraits = 25, evf = 1), fill = "#ffffff", shape = 8) +
  xlab("Number of included traits") +
  ylab("Area under receiver operating characteristic") +
  theme(legend.position = "right")

dev.off()
