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
library(ggpubr)
library(extrafont)
library(gtable)
library(grid)
library(gridExtra)

font_import()
loadfonts()

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

outputStatisticsPerTrait <- read_tsv("output/tables/outputStatisticsPerTrait.txt") %>%
  mutate(
    # trait = stringr::str_pad(trait, width = 42, side = "left", pad = " "),
    trait = ordered(trait, levels = trait),
    pearsonSquared = pearsonCorrelation^2,
    spearmanSquared = spearmanCorrelation^2,
    pearsonSquaredAdjusted = r2OnAdjustedPhenotype,
    spearmanSquaredAdjusted = spearmanCorrelationOnAdjustedPhenotype^2,
    value = case_when(traitDataType == "continuous" ~ pearsonSquared,
                      traitDataType == "ordinal" ~ spearmanCorrelation,
                      traitDataType == "binary" ~ rocCurveAuc),
    valueAdjusted = case_when(traitDataType == "continuous" ~ pearsonSquaredAdjusted,
                              traitDataType == "ordinal" ~ spearmanCorrelationOnAdjustedPhenotype,
                              traitDataType == "binary" ~ NA_real_),
    reportedValue = case_when(traitDataType == "continuous" ~ reportedR2,
                              traitDataType == "ordinal" ~ reportedR2,
                              traitDataType == "binary" ~ reportedAuc),
    reportedValueAdjusted = case_when(traitDataType == "continuous" ~ reportedR2OnAdjustedPhenotype,
                              traitDataType == "ordinal" ~ reportedR2OnAdjustedPhenotype,
                              traitDataType == "binary" ~ NA_real_),
    label = case_when(traitDataType == "continuous" ~ "squared Pearson correlation",
                      traitDataType == "ordinal" ~ "squared Spearman correlation",
                      traitDataType == "binary" ~ "area under the ROC curve")) %>%
  select(trait, traitDataType, value, valueAdjusted, reportedValue, reportedValueAdjusted, label) %>%
  pivot_longer(cols = c(value, valueAdjusted, reportedValue, reportedValueAdjusted), 
               names_to = "valueLabel", values_to = "value", values_drop_na = TRUE) %>%
  mutate(valueLabel = ordered(valueLabel, levels = c("value", "valueAdjusted", "reportedValue", "reportedValueAdjusted")))

pdf(paste0("output/figures/outputStatisticsPerTraitContinuous_", format(Sys.Date(), "%Y%m%d"), ".pdf"), 
    useDingbats = FALSE, width = 8, height = 8)

par(xpd = NA)

plotContinuous <- ggplot(
  data = outputStatisticsPerTrait %>% filter(traitDataType == "continuous" & valueLabel != "reportedValue"), 
  aes(y = value, x = trait, fill = valueLabel)) +
  geom_bar(position=position_dodge2(width = 1, preserve = "single"), stat="identity") +
  scale_fill_manual(
    values = c("#0F358C", "#007E7E", "#D19D00"),
    breaks = c("value", "valueAdjusted", "reportedValueAdjusted"),
    labels = c("UGLI, unadjusted", "UGLI, adjusted", "Literature, adjusted")) +
  ylab("Squared Pearson correlation coefficient R2") +
  xlab("Continuous traits") +
  ggtitle("Continuous traits") +
  theme(legend.position = c(0.8,0.8),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90))

plotContinuous

# plotBinary <- ggplot(
#   data = outputStatisticsPerTrait %>% filter(traitDataType == "binary"), 
#   aes(y = value, x = trait, fill = valueLabel)) +
#   geom_bar(position=position_dodge2(width = 1, preserve = "single"), stat="identity") +
#   scale_fill_manual(
#     values = c("#3980A6", "#075985"),
#     breaks = c("value", "reportedValue"),
#     labels = c("UGLI", "Literature")) +
#   ylab("Area under ROC") +
#   xlab("Binary traits") +
#   ggtitle("Binary traits") +
#   theme(legend.position = c(0.8,0.8),
#         axis.title.x = element_blank(),
#         axis.text.x = element_text(angle = 90))
# 
# plotOrdinal <- ggplot(
#   data = outputStatisticsPerTrait %>% filter(traitDataType == "ordinal"), 
#   aes(y = value, x = trait, fill = valueLabel)) +
#   geom_bar(position=position_dodge2(width = 1, preserve = "single"), stat="identity") +
#   scale_fill_manual(
#     values = c("#3980A6", "#009445", "#075985", "#39B774"),
#     breaks = c("value", "valueAdjusted", "reportedValue", "reportedValueAdjusted"),
#     labels = c("UGLI, unadjusted", "UGLI, adjusted", "Literature, unadjusted", "Literature, adjusted")) +
#   ylab("Spearman correlation") +
#   xlab("Ordinal traits") +
#   ggtitle("Ordinal traits") +
#   theme(legend.position = c(0.8,0.8),
#         axis.title.x = element_blank(),
#         axis.text.x = element_text(angle = 90))

# g2 <- ggplotGrob(plotContinuous)
# g3 <- ggplotGrob(plotOrdinal)
# g4 <- ggplotGrob(plotBinary)
# grid.arrange(g2, g3, g4, layout_matrix = rbind(c(1, 1, 1, 1, 1, 1),
#                                                  c(2, NA, 3, 3, NA, NA)),
#              top = "Title of the page")

dev.off()

