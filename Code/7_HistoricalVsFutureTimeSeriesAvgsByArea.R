# HEADER ####

# Author: Ryan Gasbarro; rygasbar@ucsc.edu
# Script Purpose: Create maps showing where projected core habitat exists 
# in National Marine Sanctuaries (NMS) in each time period
# For CCE seabird species using habitat suitability model outputs
# See Gasbarro et al. (2025) PLoS Clim. for information

# Packages ####
library(astsa, quietly=TRUE, warn.conflicts=FALSE)
library(tidyverse)
library(raster)
library(ggsci)
library(ggthemes)
library(knitr)
library(printr)
library(plyr)
library(lubridate)
library(gridExtra)
library(reshape2)
library(TTR)
library(gdata)
library(viridis)
library(RColorBrewer)
library(viridisLite)
library(wesanderson)

# Plot Theme ####
theme_Publication <- function(base_size=14, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation()
    + theme(plot.title = element_text(size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "top",
            legend.direction = "horizontal",
            legend.key.size= unit(0.5, "cm"),
            legend.spacing = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            legend.box.background = element_rect(size = 1),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_blank(),
            strip.text = element_text()
    ))
  
}

# Load Data & Wrangle ####
data.dir <- "~/SeabirdModels/Data" # ***** CHANGE TO YOUR DATA DIRECTORY *****
setwd(data.dir)

# NOTE: All time-series must be created in step 4 and combined [e.g. using rbind()] before completing this step **********************
# Alternatively, all combined time-series are provided as 'TimeSeriesAveragesByArea_AllSpecies.csv

all_cce_ts_raw <-  read.csv('TimeSeriesAveragesByArea_AllSpecies.csv')[,-1] # delete first 2 column: non-meaningful data
head(all_cce_ts_raw)

get_season <- function(date) {
  month <- month(date)
  day <- day(date)
  ifelse((month == 3 & day >= 15) | (month > 3 & month < 8) | (month == 8 & day <= 14), "Upwelling",
         ifelse((month == 8 & day >= 15) | (month > 8 & month < 11) | (month == 11 & day <= 14), "Oceanic",
                "Davidson"))
}

all_cce_ts_raw <- all_cce_ts_raw %>% mutate(season = factor(get_season(Date)))

summary_hist <- all_cce_ts_raw %>%
  filter(Date >= as.Date("1985-01-01") & Date <= as.Date("2015-12-31")) %>%
  dplyr::group_by(Species,Area, Mod, season) %>%
  dplyr::summarise(mean_prob_hist = mean(Prob, na.rm = TRUE),
            sd_prob_hist = sd(Prob, na.rm = TRUE)) %>% ungroup()

summary_fut <- all_cce_ts_raw %>%
  filter(Date >= as.Date("2070-01-01") & Date <= as.Date("2100-12-31")) %>%
  dplyr::group_by(Species,Area, Mod, season) %>%
  dplyr::summarise(mean_prob_fut = mean(Prob, na.rm = TRUE),
            sd_prob_fut = sd(Prob, na.rm = TRUE)) %>% ungroup()

summary <- merge(summary_hist, summary_fut, by = c("Species", "Area", "Mod", 'season'))
unique(summary$Area)

summary_allmod <- summary %>%
  dplyr::group_by(Species, Area, season) %>%
  dplyr::summarise(mean_prob_hist_ensemble = mean(mean_prob_hist, na.rm = TRUE),
                   sd_prob_hist_ensemble = sd(mean_prob_hist, na.rm = TRUE),
                   mean_prob_fut_ensemble = mean(mean_prob_fut, na.rm = TRUE),
                   sd_prob_fut_ensemble = sd(mean_prob_fut, na.rm = TRUE)) %>% ungroup()

head(summary_allmod)


summary$Area <- factor(summary$Area, levels = c('ShelfSlope', '200m',
                                                'CBNMS', 'GFNMS', 'MBNMS', 'CHNMS',
                                                'CoosBay', 'Brookings', 'Humboldt', 'MorroBay'))

summary(summary$Area)

summary$AreaType <- ifelse(grepl("NMS$", summary$Area), "NMS",
                           ifelse(summary$Area %in% c("MorroBay", "Humboldt", "Brookings", "CoosBay"), "WEAs",
                                  ifelse(summary$Area %in% c("200m", "ShelfSlope"), "Reference", NA)))


summary_allmod$AreaType <- ifelse(grepl("NMS$", summary_allmod$Area), "NMS",
                                  ifelse(summary_allmod$Area %in% c("MorroBay", "Humboldt", "Brookings", "CoosBay"), "WEAs",
                                         ifelse(summary_allmod$Area %in% c("200m", "ShelfSlope"), "Reference", NA)))

summary_allmod$percent_change <- (summary_allmod$mean_prob_fut_ensemble - summary_allmod$mean_prob_hist_ensemble) / summary_allmod$mean_prob_hist_ensemble * 100
summary_allmod$sd_prob_hist_ensemble
# Calculate standard errors
summary_allmod$std_err_percent_change <- sqrt(
  (100 / summary_allmod$mean_prob_hist_ensemble) ^ 2 * summary_allmod$sd_prob_fut_ensemble ^ 2 + 
    (-100 * (summary_allmod$mean_prob_fut_ensemble - summary_allmod$mean_prob_hist_ensemble) / summary_allmod$mean_prob_hist_ensemble ^ 2) ^ 2 * summary_allmod$sd_prob_hist_ensemble ^ 2
)


summary_allmod$Area <- factor(summary_allmod$Area, levels = c('ShelfSlope', '200m',
                                                              'CBNMS', 'GFNMS', 'MBNMS', 'CHNMS',
                                                              'CoosBay', 'Brookings', 'Humboldt', 'MorroBay'))
summary_allmod <- summary_allmod %>%
  group_by(Species) %>%
  mutate(
    min_prob_hist = min(mean_prob_hist_ensemble, na.rm = TRUE),
    max_prob_hist = max(mean_prob_hist_ensemble, na.rm = TRUE),
    normalized_mean_prob_hist = (mean_prob_hist_ensemble - min_prob_hist) / (max_prob_hist - min_prob_hist)
  ) %>%
  ungroup() %>%
  dplyr::select(-min_prob_hist, -max_prob_hist)
range(summary_allmod$mean_prob_hist_ensemble)

summary_allmod <- 
  summary_allmod %>%
  mutate(province = case_when(
    Area == "200m" ~ "NA",
    Area == "ShelfSlope" ~ "NA",
    Area == "CHNMS" ~ "Central",
    Area == "MBNMS" ~ "Central",
    Area == "GFNMS" ~ "Central",
    Area == "CBNMS" ~ "Central",
    Area == "MorroBay" ~ "Central",
    Area == "Humboldt" ~ "North",
    Area == "Brookings" ~ "North",
    Area == "CoosBay" ~ "North",
    TRUE ~ "NA"  # Default case
  ))



summary_allmod$province

pal4 <- wes_palette("Moonrise2")

species_labels <- c(
  "BFAL" = "Black-footed albatross",
  "CAAU" = "Cassin's auklet",
  "COMU" = "Common murre",
  "RHAU" = "Rhinoceros auklet",
  "SOSH" = "Sooty shearwater"
)

summary_allmod2$province <- factor(summary_allmod2$province, levels = c('South', 'Central', 'North', 'NA'))
summary_allmod2$season <- factor(summary_allmod2$season, levels = c('Davidson', 'Upwelling', 'Oceanic'))

# Plot ####

ggplot(summary_allmod2, aes(y = Area, x = percent_change, fill = province, alpha = normalized_mean_prob_hist)) +
  geom_errorbar(aes(xmin = percent_change - std_err_percent_change, xmax = percent_change, width = 0), show.legend = FALSE) +
  geom_bar(stat = "identity", size = 2) +
  facet_wrap(~season + Species, ncol = 5, labeller = labeller(Species = species_labels)) +
  scale_fill_manual(values = pal4) +  # Use fill for color
  labs(
    y = "Area", 
    x = "% Change", 
    fill = "Province",
    alpha = 'Hist. Suitability'
  ) +
  # guides(fill = 'none') +
  theme_Publication() + 
  geom_hline(yintercept = c(4.5, 8.5), linetype = "dashed", color = "black") +  # Add dashed lines
  scale_y_discrete(limits = rev(levels(summary_allmod2$Area))) +  
  theme(
    plot.margin = unit(c(0.5, 0.5, 0.15, 0.25), "cm"),
    legend.background = element_blank(),
    legend.box.background = element_rect(size = 0.33), 
    panel.grid.major.x = element_line(color = scales::alpha("lightgrey", 0.2)),
    strip.text = element_text(size = 0.01, color = 'white')
  )

