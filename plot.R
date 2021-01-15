#!/usr/bin/env Rscript

#############################################################
## c.a.warmerdam@umcg.nl - January 2021
## Series of functions for plotting intermediate figures
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

sampleCouplingFilePath <- "/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/data/lifelines/processed/pgs.sample-coupling-file.ugli.20210102.perm_5120samples_51mixUps.txt"


correctionModel <- lm(VALUE ~ AGE + SEX + AGE:SEX, data = completeTable)
completeTable$VALUE.CORRECTED <- resid(correctionModel)

completeTable <- completeTable %>%
  left_join(link %>% select(originalPheno = pheno, original), by=c("geno" = "original")) %>%
  mutate(
    group = case_when(geno == original ~ "null",
                      pheno == original ~ "swap",
                      TRUE ~ "alternative"),
    VALUE.CORRECTED = VALUE - predict(correctionModel, ., type = "response")
  ) %>%
  full_join(phenotypeTable %>% select(pheno, originalValue = VALUE), by = c("originalPheno" = "pheno"))

table(completeTable %>% filter(group == "swap") %>% pull(VALUE))
completeTable %>% filter(VALUE != originalValue)
  
  pdf(paste0("scatterplot_", format(Sys.time(), "%Y%m%d"), ".pdf"), 
      width=126/72, height = 117/72, useDingbats = FALSE)
  
  par(xpd = NA)
  
  # ggplot(data = completeTable, aes(x=PGS, y=VALUE.CORRECTED)) +
  #   geom_point(shape=20, color = "#009444") +
  #   geom_smooth(method=lm, se=FALSE, color="Black") +
  #   theme_classic() +
  #   theme(axis.title=element_blank(),
  #         axis.text=element_blank())
  
  ggplot(data = completeTable, aes(x=PGS, y=VALUE.CORRECTED, color=group)) +
    geom_point(shape=20) +
    scale_colour_manual(values = c("null" = "#009444", "swap" = "#D12F00")) +
    geom_smooth(data = completeTable, aes(x=PGS, y=VALUE.CORRECTED), method=lm, se=FALSE, color="Black", 
                size = (1 * ggplot2::.pt * 72.27/96)) +
    ylim(-30, NA) +
    xlim(-1, NA) +
    theme(axis.text=element_blank(),
          legend.position = "none")
  
  dev.off()

minPolygenicScore <- min(completeTable$PGS - 0.1)
maxPolygenicScore <- max(completeTable$PGS + 0.5)

pdf(paste0("hairColourSigmoid_", format(Sys.time(), "%Y%m%d"), ".pdf"),
    width=126/72, height = 117/72, useDingbats = FALSE)

par(xpd = NA)

ggplot(completeTable, aes(x=PGS, y=VALUE)) +
  geom_smooth(method = "glm", method.args = list(family="binomial"), se = FALSE,
              colour='black', size=2, fullrange=TRUE) +
  geom_histogram(data = completeTable, binwidth = 0.24, aes(x=PGS, y=..ndensity../2.4, fill=as.factor(VALUE)), position = "identity", size=0) +
  geom_point(data = completeTable %>% filter(geno != original), aes(x=PGS, y=VALUE, fill=as.factor(VALUE != originalValue)), size = (1 * ggplot2::.pt * 72.27/96), pch = 21) + 
  scale_x_continuous(limits=c(minPolygenicScore - 1, maxPolygenicScore + 1)) +
  coord_cartesian(xlim = c(minPolygenicScore, maxPolygenicScore), ylim = c(0, 1)) +
  theme(axis.text=element_blank(),
        legend.position="none")

dev.off()

plotGroupedResiduals <- function(residualsMatrix, link, filePath = NULL) {
  if (is.null(filePath)) {
    filePath <- paste0("residualsHistogram_", format(Sys.Date(), "%Y%m%d"), ".pdf")
  }

  # Convert residuals matrix to a data frame.
  residualsDataFrame <-
    as.data.frame.table(residualsMatrix, responseName = "residuals") %>%
    inner_join(link, by = c("Var1" = "pheno")) %>%
    mutate(group = case_when(Var2 == geno ~ "null",
                             TRUE ~ "alternative")) %>%
    mutate(
      diag = case_when(
        geno == Var2 ~ T,
        geno != Var2 ~ F),
      mixUp = case_when(
        geno == original & diag ~ F,
        geno != original & diag ~ T),
      groupColour = case_when(diag & mixUp ~ "red",
                              diag & !mixUp ~ "green",
                              !diag ~ "grey")
    )
  
}

pdf("resid2.pdf", 
    width=135/72, height = 117/72, useDingbats = FALSE)

par(xpd = NA)

# ggplot(scaledResidualsDataFrame, aes(x=scaledResiduals, stat(density), fill=group, size=group)) +
#   geom_histogram(bins = 32, position="identity", colour = "black") +
#   scale_fill_manual(values=c("null" = "#009444", "alternative" = "#D1D3D4")) +
#   scale_size_manual(values=c("null" = (0.5 / (ggplot2::.pt * 72.27/96)), "alternative" = 0)) +
#   geom_point(data = scaledResidualsDataFrame %>% filter(geno != original & group == "null"), aes(x=scaledResiduals, y=0), inherit.aes=F) +
#   coord_cartesian(xlim = c(-22, 22)) +
#   xlab("Residuals") +
#   theme(axis.line.y=element_blank(),axis.text.y=element_blank(),
#                 axis.ticks.y=element_blank(),
#                 axis.title.y=element_blank(),
#         legend.position = "none") +
#   stat_function(fun = dnorm, args = list(mean = mean(df$PF), sd = sd(df$PF)))

base <- ggplot(data.frame(x = c(-20, 20)), aes(x))
base + 
  stat_function(fun = dnorm, args = list(mean = mean(scaledResidualsDataFrame$scaledResiduals[scaledResidualsDataFrame$geno == scaledResidualsDataFrame$Var2]), 
                                         sd = sd(scaledResidualsDataFrame$scaledResiduals[scaledResidualsDataFrame$geno == scaledResidualsDataFrame$Var2]))) + 
  stat_function(fun = dnorm, args = list(mean = mean(scaledResidualsDataFrame$scaledResiduals[scaledResidualsDataFrame$geno != scaledResidualsDataFrame$Var2]), 
                                       sd = sd(scaledResidualsDataFrame$scaledResiduals[scaledResidualsDataFrame$geno != scaledResidualsDataFrame$Var2])))

dev.off()

pdf(paste0("residualsHistogram_", format(Sys.Date(), "%Y%m%d"), ".pdf"), 
    width=135/72, height = 117/72, useDingbats = FALSE)

par(xpd = NA)

# ggplot(residualsDataFrame, aes(x=residuals, stat(density), fill=group, size=group)) +
#   geom_histogram(bins = 16, position="identity", colour = "black") +
#   scale_fill_manual(values=c("null" = "#009444", "alternative" = "#D1D3D4")) +
#   scale_size_manual(values=c("null" = (0.5 / (ggplot2::.pt * 72.27/96)), "alternative" = 0)) +
#   geom_rug(data = residualsDataFrame %>% filter(geno != original & group == "null"), 
#              aes(x=residuals), inherit.aes=F) +
#   xlim(c(-25, 25)) +
#   # coord_cartesian(xlim = c(-22, 22)) +
#   xlab("Residuals") +
#   theme(axis.line.y=element_blank(),axis.text.y=element_blank(),
#                 axis.ticks.y=element_blank(),
#                 axis.title.y=element_blank(),
#         legend.position = "none")
  
ggplot(residualsDataFrame, aes(x=residuals, y=as.factor(groupColour), fill=groupColour)) +
  geom_violin(trim = F, draw_quantiles = c(0.25, 0.5, 0.75)) +
  scale_fill_manual(values=c("red" = "#D13B00", "green" = "#007E7E", "grey" = "#D1D3D4")) +
  geom_boxplot(size = 0.5, outlier.shape = NA) +
  xlim(c(-30, 30)) +
  # coord_cartesian(xlim = c(-22, 22)) +
  xlab("Residuals") +
  theme(axis.line.y=element_blank(),axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        legend.position = "none")

dev.off()
  
pdf("llr.pdf", 
    width=135/72, height = 117/72, useDingbats = FALSE)

par(xpd = NA)

# ggplot(scaledResidualsDataFrame, aes(x=scaledResiduals, stat(density), fill=group, size=group)) +
#   geom_histogram(bins = 32, position="identity", colour = "black") +
#   scale_fill_manual(values=c("null" = "#009444", "alternative" = "#D1D3D4")) +
#   scale_size_manual(values=c("null" = (0.5 / (ggplot2::.pt * 72.27/96)), "alternative" = 0)) +
#   geom_point(data = scaledResidualsDataFrame %>% filter(geno != original & group == "null"), aes(x=scaledResiduals, y=0), inherit.aes=F) +
#   coord_cartesian(xlim = c(-22, 22)) +
#   xlab("Residuals") +
#   theme(axis.line.y=element_blank(),axis.text.y=element_blank(),
#                 axis.ticks.y=element_blank(),
#                 axis.title.y=element_blank(),
#         legend.position = "none") +
#   stat_function(fun = dnorm, args = list(mean = mean(df$PF), sd = sd(df$PF)))



base <- ggplot(data.frame(x = c(-20, 20)), aes(x))
base + 
  stat_function(fun = function(.x) dnorm(.x, alternativeMean, alternativeSd, log = TRUE) - dnorm(.x, nullMean, nullSd, log = TRUE))

dev.off()

# at x == 20, y = 3.393156
# at x == 0, y = -0.3603929

logLikelihoodRatiosDataFrame <- 
  as.data.frame.table(logLikelihoodRatios, responseName = "logLikelihoodRatios") %>%
  inner_join(link, by = c("Var1" = "pheno"))

logLikelihoodRatiosDataFrame$group <- "alternative"
logLikelihoodRatiosDataFrame$group[logLikelihoodRatiosDataFrame$geno == logLikelihoodRatiosDataFrame$Var2] <- "null"
logLikelihoodRatiosDataFrame$group <- factor(logLikelihoodRatiosDataFrame$group, levels = c("null", "alternative"))

logLikelihoodRatiosDataFrame <- logLikelihoodRatiosDataFrame %>% 
  group_by(Var1) %>%
  mutate(scaledLlr = scale(logLikelihoodRatios)[,1])

llrSummary <- logLikelihoodRatiosDataFrame %>%
  group_by(group) %>%
  summarise(
    meanLlr = mean(logLikelihoodRatios),
    sdLlr = sd(logLikelihoodRatios)
  ) %>%
  mutate(lower = meanLlr - sdLlr,
         higher = meanLlr + sdLlr)

pdf(paste0("tTestVisual_", format(Sys.Date(), "%Y%m%d"), ".pdf"), 
    width=36*5/72, height = 117/72, useDingbats = FALSE)

par(xpd = NA)

ggplot(data = logLikelihoodRatiosDataFrame, aes(x = group, y = logLikelihoodRatios, colour = group)) +
  geom_violin(trim = TRUE, scale = "width") +
  geom_boxplot(size = 0.5, outlier.shape = NA) +
  # geom_dotplot(data = logLikelihoodRatiosDataFrame %>% filter(group == "null" & geno != original),
  #              binaxis = "y", stackdir = "center", inherit.aes = F,
  #              stackratio = 1.0, binwidth = 0.2,
  #              aes(x = group, y = logLikelihoodRatios)) +
  geom_histogram(data = logLikelihoodRatiosDataFrame %>% filter(group == "null"),
                 position="dodge", inherit.aes = F, binwidth = 0.2,
                 aes(y = logLikelihoodRatios, x=..density.., fill = as.factor(geno != original))) +
  # geom_point(data = llrSummary, aes(x = group, y = meanLlr), inherit.aes = F) +
  coord_cartesian(ylim = c(-0.3, 2.0)) + 
  # geom_linerange(data = llrSummary, aes(ymin = lower, ymax = higher, x = group), width=0.2, inherit.aes = F) +
  # geom_rug(data = logLikelihoodRatiosDataFrame %>% filter(geno != original & group == "null"), aes(y=logLikelihoodRatios), 
  #             inherit.aes=F) +
  theme(legend.position = "none")

dev.off()

# 
# dat <- tibble(group = c(rep("1", 2000), rep("2", 20)), data = c(rnorm(2000, 1, 1), rnorm(20, 4, 1)))
# 
# binwidth <- 0.3
# breaks <- sort(c(
#   seq(from = 0 - binwidth, to = min(dat$data) - binwidth, by = -binwidth), 
#   seq(from = 0, to = max(dat$data) + binwidth, by = binwidth)))
# labels <- breaks[c(2:length(breaks))] - (binwidth / 2)
# 
# dit <- dat %>%
#   mutate(
#     bin = cut(data, breaks = breaks, 
#       labels = labels, na.rm=TRUE)
#   ) %>%
#   group_by(bin, group) %>%
#   summarise(
#     count = n(),
#     logged = max(0, log10(count) + 1),
#     loggedCount = max(0, log10(count) + 1)
#   )
# 
# dit2 <- as_tibble(lapply(dit, rep, dit$loggedCount)) %>%
#   mutate(bin = as.double(as.character(bin)))

myplot <- ggplot(data = dat, aes(x = 1, y = data)) +
  geom_violin(trim = TRUE) +
  # geom_dotplot(data = dit2,
  #              method="histodot", binaxis = "y", stackdir = "center", stackgroups = TRUE, inherit.aes = F,
  #              stackratio = 1.0, binwidth = binwidth, origin = 0,
  #              aes(x = 1, y = bin, fill = group))
  geom_histogram(data = dat, position="dodge",
                 inherit.aes = F, binwidth = 0.2,
                 aes(y = data, x=..density.., fill = group))

myplot

q <- ggplot_build(myplot)
q

###################################
llrDataFrame <- as_tibble(fread("aggregatedLogLikelihoodRatiosDataFrame.tsv")) %>%
  mutate(
    diag = case_when(
      geno == Var2 ~ T,
      geno != Var2 ~ F),
    mixUp = case_when(
      geno == original & diag ~ T,
      geno != original & diag ~ F
    )) 

llrDataFrame <- llrDataFrame %>% 
  group_by(Var1) %>%
  mutate(scaledLlr = scale(logLikelihoodRatios)[,1])

densityMixUps <- density(llrDataFrame$scaledLlr[llrDataFrame$mixUp]) # set from,to equal to each other

pdf(paste0("densityRidgesWithMixup", format(Sys.Date(), "%Y%m%d"), ".pdf"), 
    width=36*5/72, height = 117/72, useDingbats = FALSE)

par(xpd = NA)

mixedUpSamples <- llrDataFrame %>% filter(diag & mixUp) %>% pull(Var1)
realSamples <- llrDataFrame %>% filter(diag & !mixUp) %>% slice_sample(n = 51) %>% pull(Var1)

ggplot(llrDataFrame %>% filter(Var1 %in% mixedUpSamples & is.finite(scaledLlr)),
       aes(y = Var1)) +
  # Add ridge lines
  geom_density_ridges_gradient(
    aes(x = scaledLlr,
        point_alpha = as.numeric(!diag),
        point_color = as.numeric(mixUp)), # Or use regular geom_point
    jittered_points = TRUE,
    position = position_points_jitter(width = 0, height = 0),
    point_shape = '|', point_size = 3, alpha = 0.7) +
  scale_colour_gradient(low = "#0F358C",
    high = "#D13B00", guide = "colourbar", aesthetics = "fill") +
  scale_point_color_continuous(high = "#0072B2", low = "#D55E00") +
  geom_vline(xintercept = c(1:2)) +
  theme(legend.position = "none")
  
ggplot(llrDataFrame %>% filter(Var1 %in% realSamples & is.finite(scaledLlr)),
       aes(y = Var1)) +
  # Add ridge lines
  geom_density_ridges(
    aes(x = scaledLlr,
        point_alpha = as.numeric(!diag),
        point_color = as.numeric(mixUp)),
    jittered_points = TRUE,
    position = position_points_jitter(width = 0, height = 0),
    point_shape = '|', point_size = 3, alpha = 0.7) +
  scale_point_color_continuous(high = "#0072B2", low = "#D55E00") +
  geom_vline(xintercept = c(1:2)) +
  theme(legend.position = "none")

dev.off()

mat <- fread("/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/output/sample-swap-prediction/20200811.20201008.efi-discretization.30/aggregatedLogLikelihoodRatiosMatrix.tsv")
aggregatedLlrMatrix <- as.matrix(mat[,2:ncol(mat)])

rownames(aggregatedLlrMatrix) <- mat$V1

lrProducts <- 
  as.data.frame.table(aggregatedLlrMatrix, responseName = "logLikelihoodRatios") %>%
  inner_join(link, by = c("Var1" = "pheno"))

lrProducts$group <- "alternative"
lrProducts$group[lrProducts$geno == lrProducts$Var2] <- "null"
lrProducts$group <- factor(lrProducts$group, c("alternative", "null"))

phenoSamples <- sample(unique(lrProducts[lrProducts$geno == lrProducts$Var2, "Var1"]), size = 26)
phenoSamples <- unique(lrProducts$Var1[lrProducts$geno != lrProducts$original])

mat <- fread("/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/output/sample-swap-prediction/20200811.20201012.efi-discretization.30/aggregatedLogLikelihoodRatiosMatrix.tsv")

sampledLrProducts <- lrProducts %>%
  #inner_join(numberOfTraits, by=c("Var1", "Var2")) %>%
  mutate(colourGroup = case_when(
    Var2 == geno ~ 1,
    Var2 == original ~ 2,
    TRUE ~ 3)) %>%
  filter(Var1 %in% phenoSamples)

# sampledLrProducts <- sampledLrProducts %>% 
#   group_by(Var1) %>%
#   mutate(scaledLlr = scale(logLikelihoodRatios))

cols = c("2" = "blue", "1" = "red", "3" = "black")

ggplot(sampledLrProducts %>% filter(is.finite(scaledLlr)),
       aes(y = Var1)) +

  # Add ridge lines
  geom_density_ridges(
    aes(x = scaledLlr,
        point_alpha = as.numeric(colourGroup != 3),
        point_color = colourGroup),
    jittered_points = TRUE,
    position = position_points_jitter(width = 0, height = 0),
    point_shape = '|', point_size = 3, alpha = 0.7) +
  scale_point_color_continuous(high = "#0072B2", low = "#D55E00") +
  #scale_discrete_manual(values = cols, aesthetics = "point_color") +
  #scale_discrete_manual("point_color", values = cols) +

  # Add annotation with number of traits, and the number of filtered likelihood ratios
  # geom_text(
  #   data=sampledLrProducts %>% group_by(Var1) %>%
  #     summarise(numberOfTraits = median(numberOfTraits),
  #               numberOfFiltered = sum(!is.finite(logLikelihoodRatios))),
  #   position=position_nudge(y=0.64), colour="red", size=3.5,
  #   hjust = "inward", x = 0,
  #   aes(label = sprintf("n traits: %d, n filtered: %d", numberOfTraits, numberOfFiltered))) +

  # Set theme
  theme_minimal() +
  theme(legend.position = "none")

ggsave("ridges_mixups26_scaled_height.png", width=8, height=20)