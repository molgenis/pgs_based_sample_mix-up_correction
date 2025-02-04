#!/usr/bin/env Rscript

#############################################################
## c.a.warmerdam@umcg.nl - January 2021
## Plot predictive power for identifying mix-ups
#############################################################

##############################
# Load libraries
##############################
library(tidyverse)
library(argparse)
library(data.table)
library(pROC)

##############################
# Define argument parser
##############################

parser <- ArgumentParser(description='')
parser$add_argument('--dir', required = T,
                    help='path from where to read sample swap prediction results.')
parser$add_argument('--phenotypes-file', required = T,
                    help='path to a tab-delimited file holding all processed phenotype data.')

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

# Parse arguments.
args <- parser$parse_args(commandArgs(trailingOnly = TRUE))
dir <- args$dir
phenotypesFilePath <- args$phenotypes_file

# Read aggregated likelihood ratio dataframe as well as aggregated number of traits in a matrix.
llrDataFrame <- as_tibble(fread(file.path(dir, "aggregatedLogLikelihoodRatiosDataFrame.tsv")))
numberOfTraits <- as.matrix(fread(file.path(dir, "aggregatedNumberOfTraitsMatrix.tsv")), rownames = 1)

llrDataFrame$group <- "alternative"
llrDataFrame$group[llrDataFrame$original == llrDataFrame$Var2] <- "null"
llrDataFrame$group <- ordered(llrDataFrame$group, levels = c("null", "alternative"))

llrDataFrame <- llrDataFrame %>% 
  group_by(Var1) %>%
  mutate(scaledLlr = scale(logLikelihoodRatios)[,1])

# Confine ourselves to the diagonal
permutationTestDataFrame <- llrDataFrame %>%
  filter(geno == Var2)

# Assign the number of traits to the permutationTestDataFrame
permutationTestDataFrame$numberOfTraits <- diag(numberOfTraits)

reportedSexTable <- fread(phenotypesFilePath, header=T, quote="", sep="\t") %>%
  rename_all(recode, "UGLI_ID" = "ID") %>%
  mutate(SEX = factor(SEX, levels = c("Female", "Male"))) %>%
  group_by(ID) %>%
  filter(!any(AGE < 18)) %>%
  select(ID, SEX) %>%
  distinct()

stopifnot(nrow(reportedSexTable) == length(unique(reportedSexTable$ID)))

# Assuming a sex check has already been performed on the data, the reported sex should
# already correspond to the inferred sex.
sexCheckData <- reportedSexTable %>% rename_all(recode, "ID" = "UGLI_ID")
levels(sexCheckData$SEX) <- list(F = "Female", M = "Male")
sexCheckData$Reported_sex <- sexCheckData$SEX
sexCheckData$Inferred_sex <- sexCheckData$SEX

# Check correspondence.
all(sexCheckData$Reported_sex == sexCheckData$Inferred_sex)

# Assign the sex concordance data
permutationTestDataFrame <- llrDataFrame %>%
  filter(geno == Var2) %>%
  inner_join(sexCheckData %>% select(UGLI_ID, Reported_sex), by = c("Var1" = "UGLI_ID")) %>%
  inner_join(sexCheckData %>% select(UGLI_ID, Inferred_sex), by = c("geno" = "UGLI_ID")) %>%
  mutate(sexCheck = case_when(Reported_sex == Inferred_sex ~ 0L,
                              Reported_sex != Inferred_sex ~ 1L)) %>%
  
  # The when the sexCheck identifies a mix-up, we set the predicted value to infinity.
  mutate(combined = case_when(sexCheck == 1 ~ Inf,
                              TRUE ~ scaledLlr))

# Calculate ROC for PGS based method only
rocPgs <- roc(permutationTestDataFrame$group ~ permutationTestDataFrame$scaledLlr)
# Calculate ROC for sex check only
rocSexCheck <- roc(permutationTestDataFrame$group ~ permutationTestDataFrame$sexCheck)
# Calculate ROC for combined model
rocCombined <- roc(permutationTestDataFrame$group ~ permutationTestDataFrame$combined)


roclist <- list("PGS based method" = rocPgs,
                "Sex-check" = rocSexCheck,
                "Combined" = rocCombined)

# Loop through ROC thresholds for reporting data.
fdrThresholds <- c(0.0005, 0.001, 0.005, 0.01, 0.05, 0.25, 0.5, 0.9)

# For each of the ROCs, generate confusion matrix values.
specificitiesPgs <- as_tibble(t(sapply(fdrThresholds, function(falsePositiveRate) {
  rocResults <- coords(rocPgs, 1 - falsePositiveRate, input = "specificity", 
                       ret=c("threshold", "tpr", "fpr", "tnr", "fnr"), transpose = TRUE)
  rocResults["threshold"] <- rocPgs$thresholds[which.min(abs(rocPgs$specificities - (1 - falsePositiveRate)))]
  rocResults["deviance"] <- min(abs(rocPgs$specificities - (1 - falsePositiveRate)))
  return(rocResults)
})))

specificitiesSexCheck <- as_tibble(t(sapply(fdrThresholds, function(falsePositiveRate) {
  rocResults <- coords(rocSexCheck, 1 - falsePositiveRate, input = "specificity", 
                       ret=c("threshold", "tpr", "fpr", "tnr", "fnr"), transpose = TRUE)
  rocResults["threshold"] <- rocSexCheck$thresholds[which.min(abs(rocSexCheck$specificities - (1 - falsePositiveRate)))]
  rocResults["deviance"] <- min(abs(rocSexCheck$specificities - (1 - falsePositiveRate)))
  return(rocResults)
})))

specificitiesCombined <- as_tibble(t(sapply(fdrThresholds, function(falsePositiveRate) {
  rocResults <- coords(rocCombined, 1 - falsePositiveRate, input = "specificity", 
                       ret=c("threshold", "tpr", "fpr", "tnr", "fnr"), transpose = TRUE)
  rocResults["threshold"] <- rocCombined$thresholds[which.min(abs(rocCombined$specificities - (1 - falsePositiveRate)))]
  rocResults["deviance"] <- min(abs(rocCombined$specificities - (1 - falsePositiveRate)))
  
  return(rocResults)
})))

# Generate confusion matrix values for the perfect chance diagonal
specificitiesCoinFlip <- tibble(threshold = NA_real_, fpr = fdrThresholds, tnr = 1 - fdrThresholds, tpr = fdrThresholds, fnr = 1 - fdrThresholds)

# merge confusion tabels.
specificitiesBase <- bind_rows(specificitiesPgs, specificitiesSexCheck, specificitiesCombined, specificitiesCoinFlip, .id = "rocType") %>%
  mutate(falsePositiveRate = round(fpr, digits = 8)) %>%
  mutate(rocType = ordered(case_when(rocType == 1 ~ "pgsOnly", 
                                     rocType == 2 ~ "sexCheck",
                                     rocType == 3 ~ "combined",
                                     rocType == 4 ~ "coinFlip"), c("pgsOnly", "sexCheck", "combined", "coinFlip")))

# Get increase for all thresholds.
specificitiesIncreasesPgsOnly <- specificitiesBase%>%
  mutate(fpr = as.character(fpr)) %>% 
  filter(rocType == "pgsOnly") %>%
  inner_join(specificitiesBase %>% filter(rocType == "coinFlip"), by = c("fpr" = "fpr"))

specificitiesIncreasesCombined <- specificitiesBase%>%
  mutate(fpr = as.character(fpr)) %>% 
  filter(rocType == "combined") %>%
  inner_join(specificitiesBase %>% filter(rocType == "sexCheck"), by = "fpr")

specificitiesIncreases <- bind_rows(specificitiesIncreasesPgsOnly, specificitiesIncreasesCombined) %>%
  mutate(tprIncrease = (tpr.x - tpr.y) / tpr.y * 100,
         tprIncreasePp = tpr.x - tpr.y)

print(specificitiesIncreases, width = Inf)

# Plot the ROC curves.
g <- ggroc(roclist, aes = "linetype", legacy.axes = TRUE) +
  geom_abline() +
  ggtitle("ROC curve") +
  labs(x = "False positive rate",
       y = "True positive rate",
       linetype = "Different legend title") +
  geom_point(data = specificitiesBase, aes(x = fpr, y = tpr, colour = rocType), inherit.aes = F) +
  theme(legend.position = c(0.8,0.2))

# Loop over fdr thresholds to plot these.
for (fdrThreshold in fdrThresholds) {
  
  g <- g + 
    geom_vline(xintercept = fdrThresholds)
}

# Plot the ROCs
pdf(file.path(dir, paste0("combinedRoc_", format(Sys.Date(), "%Y%m%d"), ".pdf")), 
    useDingbats = FALSE, width = 4, height = 4)

par(xpd = NA)

g

dev.off()

# Plot the confusion matrices.
specificities <- specificitiesBase %>%
  mutate(falsePositiveRate = round(fpr, digits = 8)) %>%
  pivot_longer(cols = c("tpr", "fpr", "tnr", "fnr"), names_to = "statistic", values_to = "value") %>%
  mutate(trueClass = case_when(statistic %in% c("tpr", "fnr") ~ "actual mix-up",
                               statistic %in% c("fpr", "tnr") ~ "not actual mix-up"),
         predictedClass = case_when(statistic %in% c("tpr", "fpr") ~ "predicted mix-up",
                                    statistic %in% c("fnr", "tnr") ~ "not predicted mix-up"))

# Plot the confusion matrices
pdf(file.path(dir, paste0("confusionMatrices_", format(Sys.Date(), "%Y%m%d"), ".pdf")), 
    useDingbats = FALSE, width = 8, height = 11.5)

par(xpd = NA)

ggplot(data =  specificities, mapping = aes(x = trueClass, y = predictedClass)) +
  geom_tile(aes(fill = value), colour = "white", size = 0, width = 1, height = 1) +
  geom_text(aes(label = sprintf("%1.2f", value)), vjust = "center", size = 0.5 * 72.27/96 * 8) +
  scale_fill_gradient(low = "white", high = "#949494") +
  coord_equal() +
  labs(title = "Confusion matrices") +
  scale_y_discrete(expand = c(0, 0), drop = T) +
  scale_x_discrete(expand = c(0, 0), drop = T) +
  facet_grid(falsePositiveRate ~ rocType) +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_text(size = 8),
        axis.text.x = element_text(angle = 90),
        legend.position="bottom")

dev.off()

# Plot the distributions
pdf(file.path(dir, paste0("outputDistributions_", format(Sys.Date(), "%Y%m%d"), ".pdf")), 
    width=36*4/72, height = 117/72, useDingbats = FALSE)

par(xpd = NA)

ggplot(data = permutationTestDataFrame, aes(x = scaledLlr, fill = group)) +
  geom_density() +
  theme(legend.position = "none")

dev.off()
