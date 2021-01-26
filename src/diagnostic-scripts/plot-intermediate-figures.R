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
library(pROC)

# Attempts to set the default ggplot theme to classic, and update the
# axes to have a width of 1 in illustrator points.
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

# Function that plots the corrected relationship between a phenotype and PGS.
plotRelationshipContinuousTrait <- function(completeTable, link) {
  
  # Correct the phenotype
  correctionModel <- lm(VALUE ~ AGE + SEX + AGE:SEX, data = completeTable)
  completeTable$VALUE.CORRECTED <- resid(correctionModel)

  # Map point groups
  completeTable <- completeTable %>%
    left_join(link %>% select(originalPheno = pheno, geno, original), by=c("geno" = "geno")) %>%
    mutate(
      group = case_when(geno == original ~ "null",
                        pheno == original ~ "swap",
                        TRUE ~ "alternative"),
      VALUE.CORRECTED = VALUE - predict(correctionModel, ., type = "response")
    )
  
  pdf(paste0("scatterplot_", format(Sys.time(), "%Y%m%d"), ".pdf"), 
      width=126/72, height = 117/72, useDingbats = FALSE)
  
  par(xpd = NA)
  
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
  
}

# Function that plots relationship for binary traits.
plotRelationshipBinaryTrait <- function(completeTable, filePath = NULL) {
  if (is.null(filePath)) {
    filePath <- paste0("hairColourSigmoid_", format(Sys.time(), "%Y%m%d"), ".pdf")
  }
  
  minPolygenicScore <- min(completeTable$PGS - 0.1)
  maxPolygenicScore <- max(completeTable$PGS + 0.5)
  
  pdf(filePath,
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
}

# Function obtains residuals data frame from matrix and link tibble.
getResidualsDataFrame <- function(residualsMatrix, link, filePath = NULL) {
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
  
  return(residualsDataFrame)
  
}

# Function that plots two Gaussian functions given the distribution of residuals.
plotGaussianFunctions <- function(residualsDataFrame, filePath = NULL) {
  if (is.null(filePath)) {
    filePath <- paste0("residualGaussians_", format(Sys.Date(), "%Y%m%d"), ".pdf")
  }
  pdf(filePath, width=135/72, height = 117/72, useDingbats = FALSE)
  
  par(xpd = NA)

  base <- ggplot(data.frame(x = c(-20, 20)), aes(x))
  base + 
    stat_function(fun = dnorm, args = list(mean = mean(residualsDataFrame$scaledResiduals[residualsDataFrame$geno == residualsDataFrame$Var2]), 
                                           sd = sd(residualsDataFrame$scaledResiduals[residualsDataFrame$geno == residualsDataFrame$Var2]))) + 
    stat_function(fun = dnorm, args = list(mean = mean(residualsDataFrame$scaledResiduals[residualsDataFrame$geno != residualsDataFrame$Var2]), 
                                           sd = sd(residualsDataFrame$scaledResiduals[residualsDataFrame$geno != residualsDataFrame$Var2])))
  
  dev.off()
}

# Function that plots the distribution of residuals for permuted, 
# supposedly correct, and induced mix-ups.
plotResiduals <- function(residualsDataFrame, filePath = NULL) {
  if (is.null(filePath)) {
    filePath <- paste0("residualsHistogram+Violins_", format(Sys.Date(), "%Y%m%d"), ".pdf")
  }
  
  pdf(filePath, width=135/72, height = 117/72, useDingbats = FALSE)
  
  par(xpd = NA)
  
  ggplot(residualsDataFrame, aes(x=residuals, stat(density), fill=group, size=group)) +
    geom_histogram(bins = 16, position="identity", colour = "black") +
    scale_fill_manual(values=c("null" = "#009444", "alternative" = "#D1D3D4")) +
    scale_size_manual(values=c("null" = (0.5 / (ggplot2::.pt * 72.27/96)), "alternative" = 0)) +
    geom_rug(data = residualsDataFrame %>% filter(geno != original & group == "null"),
               aes(x=residuals), inherit.aes=F) +
    xlim(c(-25, 25)) +
    # coord_cartesian(xlim = c(-22, 22)) +
    xlab("Residuals") +
    theme(axis.line.y=element_blank(),axis.text.y=element_blank(),
                  axis.ticks.y=element_blank(),
                  axis.title.y=element_blank(),
          legend.position = "none")
    
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
}

# Plots the likelihood ratio curve given two fitted gaussian functions.
plotLogLikelihoodRatioCurve <- function(gaussianParameters, filePath = NULL) {
  if (is.null(filePath)) {
    filePath <- paste0("logLikelihoodRatioCurve_", format(Sys.Date(), "%Y%m%d"), ".pdf")
  }
  
  # Calculate the mean of both groups of residuals.
  nullMean <- gaussianParameters$null$mean
  alternativeMean <- gaussianParameters$alternative$mean
  
  # Calculate the standard deviation in both groups of residuals.
  nullSd <- gaussianParameters$null$sd
  alternativeSd <- gaussianParameters$alternative$sd
  
  pdf("llr.pdf", 
      width=135/72, height = 117/72, useDingbats = FALSE)
  
  par(xpd = NA)
  
  base <- ggplot(data.frame(x = c(-20, 20)), aes(x))
  base + 
    stat_function(fun = function(.x) dnorm(.x, alternativeMean, alternativeSd, log = TRUE) - dnorm(.x, nullMean, nullSd, log = TRUE))
  
  dev.off()
}

# Plot two distributions and a histogram comparing the provided and permuted samples,
# as well as the supposedly correct and swapped samples.
plotTtestVisual <- function(logLikelihoodRatios, link, filePath = NULL) {
  if (is.null(filePath)) {
    filePath <- paste0("tTestVisual_", format(Sys.Date(), "%Y%m%d"), ".pdf")
  }
  
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
  
  pdf(filePath, width=36*5/72, height = 117/72, useDingbats = FALSE)
  par(xpd = NA)
  
  ggplot(data = logLikelihoodRatiosDataFrame, aes(x = group, y = logLikelihoodRatios, colour = group)) +
    geom_violin(trim = TRUE, scale = "width") +
    geom_boxplot(size = 0.5, outlier.shape = NA) +
    geom_histogram(data = logLikelihoodRatiosDataFrame %>% filter(group == "null"),
                   position="dodge", inherit.aes = F, binwidth = 0.2,
                   aes(y = logLikelihoodRatios, x=..density.., fill = as.factor(geno != original))) +
    coord_cartesian(ylim = c(-0.3, 2.0)) + 
    theme(legend.position = "none")
  
  dev.off()
}