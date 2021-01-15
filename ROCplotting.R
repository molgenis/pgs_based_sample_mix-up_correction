#!/usr/bin/env Rscript

#############################################################
## c.a.warmerdam@umcg.nl - January 2021
## Plot predictive power of PGSes for traits
#############################################################

##############################
# Load libraries
##############################
library(tidyverse)
library(data.table)
library(RColorBrewer)
library(pROC)

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

llrDataFrame <- as_tibble(fread("aggregatedLogLikelihoodRatiosDataFrame.tsv"))
numberOfTraits <- as.matrix(fread("aggregatedNumberOfTraitsMatrix.tsv"), rownames = 1)

llrDataFrame$group <- "alternative"
llrDataFrame$group[llrDataFrame$original == llrDataFrame$Var2] <- "null"
llrDataFrame$group <- ordered(llrDataFrame$group, levels = c("null", "alternative"))

llrDataFrame <- llrDataFrame %>% 
  group_by(Var1) %>%
  mutate(scaledLlr = scale(logLikelihoodRatios)[,1])

# Confine ourselves to the diagonal
permutationTestDataFrame <- llrDataFrame %>%
  filter(geno == Var2)

permutationTestDataFrame$numberOfTraits <- diag(numberOfTraits)

phenotypesFilePath <- "/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/data/lifelines/processed/pgs.phenotypes_20201215.ugli.dat"
reportedSexTable <- fread(phenotypesFilePath, header=T, quote="", sep="\t") %>%
  rename_all(recode, "UGLI_ID" = "ID") %>%
  mutate(SEX = factor(SEX, levels = c("Female", "Male"))) %>%
  group_by(ID) %>%
  filter(!any(AGE < 18)) %>%
  select(ID, SEX) %>%
  distinct()

stopifnot(nrow(reportedSexTable) == length(unique(reportedSexTable$ID)))

ugliQc <- read_tsv("/groups/umcg-lifelines/tmp01/releases/gsa_genotypes/v1/Logs/ugli_qc_release1_samples.csv")
gsaLinkage <- read_tsv("/groups/umcg-lifelines/tmp01/releases/gsa_linkage_files/v1/gsa_linkage_file.dat")
link <- read_tsv("/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/data/lifelines/processed/pgs.sample-coupling-file.ugli.20201229-other.perm_16409samples_8204mixUps.txt")

sexCheckData <- ugliQc %>%
  inner_join(gsaLinkage, by = c("Sample_ID" = "genotyping_name")) %>%
  select(UGLI_ID, Inferred_sex) %>%
  inner_join(reportedSexTable, by = c("UGLI_ID" = "ID"))
levels(sexCheckData$SEX) <- list(F = "Female", M = "Male")
sexCheckData$Reported_sex <- sexCheckData$SEX

all(sexCheckData$SEX == sexCheckData$Inferred_sex)

permutationTestDataFrame <- llrDataFrame %>%
  filter(geno == Var2) %>%
  inner_join(sexCheckData %>% select(UGLI_ID, Reported_sex), by = c("Var1" = "UGLI_ID")) %>%
  inner_join(sexCheckData %>% select(UGLI_ID, Inferred_sex), by = c("geno" = "UGLI_ID")) %>%
  mutate(sexCheck = case_when(Reported_sex == Inferred_sex ~ 0L,
                              Reported_sex != Inferred_sex ~ 1L)) %>%
  mutate(combined = case_when(sexCheck == 1 ~ Inf,
                              TRUE ~ scaledLlr))

rocPgs <- roc(permutationTestDataFrame$group ~ permutationTestDataFrame$scaledLlr)
rocSexCheck <- roc(permutationTestDataFrame$group ~ permutationTestDataFrame$sexCheck)
rocCombined <- roc(permutationTestDataFrame$group ~ permutationTestDataFrame$combined)
roclist <- list("PGS based method" = rocPgs,
                "Sex-check" = rocSexCheck,
                "Combined" = rocCombined)

fdrThresholds <- c(0.01, 0.05, 0.25, 0.5, 0.9)

specificitiesPgs <- as_tibble(t(sapply(fdrThresholds, function(fpr) {
  coords(rocPgs, 1 - fpr, input = "specificity", ret=c("threshold", "tpr", "fpr", "tnr", "fnr"), transpose = TRUE)
})))

specificitiesSexCheck <- as_tibble(t(sapply(c(0, fdrThresholds), function(fpr) {
  coords(rocSexCheck, 1 - fpr, input = "specificity", ret=c("threshold", "tpr", "fpr", "tnr", "fnr"), transpose = TRUE)
})))

specificitiesCombined <- as_tibble(t(sapply(fdrThresholds, function(fpr) {
  coords(rocCombined, 1 - fpr, input = "specificity", ret=c("threshold", "tpr", "fpr", "tnr", "fnr"), transpose = TRUE)
})))

specificitiesCoinFlip <- tibble(threshold = NA_real_, fpr = fdrThresholds, tnr = 1 - fdrThresholds, tpr = fdrThresholds, fnr = 1 - fdrThresholds)

specificitiesBase <- bind_rows(specificitiesPgs, specificitiesSexCheck, specificitiesCombined, specificitiesCoinFlip, .id = "rocType") %>%
  mutate(rocType = ordered(case_when(rocType == 1 ~ "pgsOnly", 
                                     rocType == 2 ~ "sexCheck",
                                     rocType == 3 ~ "combined",
                                     rocType == 4 ~ "coinFlip"), c("pgsOnly", "sexCheck", "combined", "coinFlip"))) %>%
  mutate(fpr = as.character(fpr))

specificitiesIncreasesPgsOnly <- specificitiesBase %>% 
  filter(rocType == "pgsOnly") %>%
  inner_join(specificitiesBase %>% filter(rocType == "coinFlip"), by = c("fpr" = "fpr"))

specificitiesIncreasesCombined <- specificitiesBase %>% 
  filter(rocType == "combined") %>%
  inner_join(specificitiesBase %>% filter(rocType == "sexCheck"), by = "fpr")

specificitiesIncreases <- bind_rows(specificitiesIncreasesPgsOnly, specificitiesIncreasesCombined) %>%
  mutate(tprIncrease = (tpr.x - tpr.y) / tpr.y * 100,
         tprIncreasePp = tpr.x - tpr.y)

print(specificitiesIncreases, width = Inf)
  
  # mutate(rocType = ordered(case_when(rocType == 1 ~ "pgsOnly", 
  #                                    rocType == 2 | rocType == 3 ~ "combined"), c("pgsOnly", "combined")),
  #        falsePositiveRate = fpr)

g <- ggroc(roclist, aes = "linetype", legacy.axes = TRUE) +
  geom_abline() +
  ggtitle("ROC curve") +
  labs(x = "False positive rate",
       y = "True positive rate",
       linetype = "Different legend title") +
  geom_point(data = specificitiesBase, aes(x = fpr, y = tpr, colour = rocType), inherit.aes = F) +
  theme(legend.position = c(0.8,0.2))

for (fdrThreshold in fdrThresholds) {
  
  g <- g + 
    geom_vline(xintercept = fdrThresholds)
}

pdf(paste0("combinedRoc_", format(Sys.Date(), "%Y%m%d"), ".pdf"), 
    useDingbats = FALSE, width = 8, height = 4)

par(xpd = NA)

g

dev.off()

specificities <- specificitiesBase %>%
  pivot_longer(cols = c("tpr", "fpr", "tnr", "fnr"), names_to = "statistic", values_to = "value") %>%
  mutate(trueClass = case_when(statistic %in% c("tpr", "fnr") ~ "actual mix-up",
                               statistic %in% c("fpr", "tnr") ~ "not actual mix-up"),
         predictedClass = case_when(statistic %in% c("tpr", "fpr") ~ "predicted mix-up",
                                    statistic %in% c("fnr", "tnr") ~ "not predicted mix-up"))

pdf(paste0("confusionMatrices_", format(Sys.Date(), "%Y%m%d"), ".pdf"), 
    useDingbats = FALSE, width = 4, height = 8)

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


pdf(paste0("outputDistributions", format(Sys.Date(), "%Y%m%d"), ".pdf"), 
    width=36*4/72, height = 117/72, useDingbats = FALSE)

par(xpd = NA)

ggplot(data = permutationTestDataFrame, aes(x = scaledLlr, fill = group)) +
  # geom_dotplot(data = logLikelihoodRatiosDataFrame %>% filter(group == "null" & geno != original),
  #              binaxis = "y", stackdir = "center", inherit.aes = F,
  #              stackratio = 1.0, binwidth = 0.2,
  #              aes(x = group, y = logLikelihoodRatios)) +
  geom_density() +
  # geom_point(data = llrSummary, aes(x = group, y = meanLlr), inherit.aes = F) +
  # coord_cartesian(xlim = c(-0.3, 2.0)) + 
  # geom_linerange(data = llrSummary, aes(ymin = lower, ymax = higher, x = group), width=0.2, inherit.aes = F) +
  # geom_rug(data = logLikelihoodRatiosDataFrame %>% filter(geno != original & group == "null"), aes(y=logLikelihoodRatios), 
  #             inherit.aes=F) +
  theme(legend.position = "none")

dev.off()
