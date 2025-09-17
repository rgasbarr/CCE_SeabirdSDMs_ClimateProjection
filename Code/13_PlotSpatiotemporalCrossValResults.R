# HEADER ####
# Author: Ryan Gasbarro; rygasbar@ucsc.edu
# Purpose: Plot spatiotemporal cross-validation (see Step 12a & 12b) Results
# See Gasbarro et al. (2025) PLoS Climate for more information

library(tidyverse)
library(grid)
library(gridExtra)

# Plot cross-validation outputs by species ####
data.dir <- "~/SeabirdModels/Data" # * CHANGE TO YOUR DATA DIRECTORY
setwd(data.dir)

bfal_loyo <- read.csv('BFAL_EvalLOYO.csv')
caau_loyo <- read.csv('CAAU_EvalLOYO.csv')
comu_loyo <- read.csv('COMU_EvalLOYO.csv')
rhau_loyo <- read.csv('RHAU_EvalLOYO.csv')
sosh_loyo <- read.csv('SOSH_EvalLOYO.csv')

bfal_loyo$Species <- 'BFAL'
caau_loyo$Species <- 'CAAU'
comu_loyo$Species <- 'COMU'
rhau_loyo$Species <- 'RHAU'
sosh_loyo$Species <- 'SOSH'


loyo <- rbind(bfal_loyo, caau_loyo, comu_loyo, rhau_loyo, sosh_loyo)
loyo <- loyo[,-1]

rm(bfal_loyo, caau_loyo, comu_loyo, rhau_loyo, sosh_loyo)

bfal_lopo <- read.csv('BFAL_Eval_LOPO.csv')
caau_lopo <- read.csv('CAAU_Eval_LOPO.csv')
comu_lopo <- read.csv('COMU_EvalLOPO.csv')
rhau_lopo <- read.csv('RHAU_EvalLOPO.csv')
sosh_lopo <- read.csv('SOSH_EvalLOPO.csv')

bfal_lopo$Species <- 'BFAL'
caau_lopo$Species <- 'CAAU'
comu_lopo$Species <- 'COMU'
rhau_lopo$Species <- 'RHAU'
sosh_lopo$Species <- 'SOSH'

lopo <- rbind(bfal_lopo, caau_lopo, comu_lopo, rhau_lopo, sosh_lopo)
lopo <- lopo[,-1]

rm(bfal_lopo, caau_lopo, comu_lopo, rhau_lopo, sosh_lopo)

my.theme <- theme(
  panel.background = element_blank(),
  panel.border = element_rect(colour="black", fill = NA),
  panel.grid = element_blank(),
  axis.line = element_line("black"),
  text = element_text(size = 10),
  axis.text = element_text(size = 9, colour = "black"),
  axis.title = element_text(size = 11, colour = "black"),
  legend.title = element_text(size = 11),
  legend.text = element_text(size = 10),
  legend.key = element_blank(),
  legend.background = element_blank(),
  legend.position = 'top')

loyo_auc <- ggplot(loyo, aes(x = k, y = AUC, col = Species, group = Species)) +
  geom_point() + geom_line() + facet_wrap(~Species, ncol = 5) + my.theme + labs(x = NULL) + guides(col = 'none')

loyo_tss <- ggplot(loyo, aes(x = k, y = TSS, col = Species, group = Species)) +
  geom_point() + geom_line() + facet_wrap(~Species, ncol = 5) + my.theme + labs(x = NULL) + guides(col = 'none')

lopo$k <- ordered(lopo$k, levels = c("SOUTH", "CENTRAL", "NORTH"))

lopo_auc <- ggplot(lopo, aes(x = k, y = AUC, col = Species, group = Species)) +
  geom_point() + geom_line() + facet_wrap(~Species, ncol = 5) + my.theme + labs(x = NULL) + guides(col = 'none')

lopo_tss <- ggplot(lopo, aes(x = k, y = TSS, col = Species, group = Species)) +
  geom_point() + geom_line() + facet_wrap(~Species, ncol = 5) + my.theme + labs(x = NULL) + guides(col = 'none')

grid.arrange(loyo_auc, loyo_tss, lopo_auc, lopo_tss, ncol = 1)

































