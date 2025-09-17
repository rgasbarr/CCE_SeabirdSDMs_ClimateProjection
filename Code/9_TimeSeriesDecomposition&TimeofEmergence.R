# HEADER ####
# Author: Ryan Gasbarro; rygasbar@ucsc.edu
# Script Purpose: Example script showing time-series decomposition and time-of-emergence calculations for
# CCE seabird species using habitat suitability model outputs
# See Gasbarro et al. (2025) PLoS Clim. for information

# Load Packages ####
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
library(gdata)
library(viridis)
library(viridisLite)

# PLOT THEMES ####
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
            panel.grid.major.x = element_line(colour = 'lightgrey', linewidth = 0.1),
            panel.grid.minor.x = element_line(colour = 'lightgrey', linewidth = 0.1),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "top",
            legend.direction = "horizontal",
            legend.key.size= unit(1, "cm"),
            legend.spacing = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            legend.box.background = element_rect(size = 1),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_blank(),
            strip.text = element_text()
    ))
  
}
scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}
scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

# _________________________ ####
# Time-series decomposition ####
# **** NOTE **** This example shows the time-series decomposition for the shelf (200m) and shelf-slope study domain. 
# ************** This process must be repeated for each Species * Area * Earth-system model for which time-series were created in Step 4
# ************** These combined time-series are provided in the file: 'TimeSeriesRawByCCEArea_AllSpecies.csv'
# IPSL 
#  bfal 200m (on-shelf #
ipsl.dir = "bfal/dailyP_IPSL/"  # change to the folder where IPSL daily SDM projections are stored
setwd(ipsl.dir)

bfal_ipsl_200m_ts_raw <- read.csv("bfal43_MeanOccProb_TimeSeries_1980-2100_200m.csv")
head(bfal_ipsl_200m_ts_raw)
bfal_ipsl_200m_ts_raw$Date <- as.Date(bfal_ipsl_200m_ts_raw$Date)
bfal_ipsl_200m_ts_raw <- bfal_ipsl_200m_ts_raw[!(month(bfal_ipsl_200m_ts_raw$Date) == 2 & day(bfal_ipsl_200m_ts_raw$Date) == 29), ] # remove leap day obs
bfal_ipsl_200m_ts_raw <- bfal_ipsl_200m_ts_raw[,-1]
bfal_ipsl_200m_ts_raw$Mod <- 'IPSL'
bfal_ipsl_200m_ts_raw$Area <- '200m'
bfal_ipsl_200m_ts_raw$AreaType <- 'CCE'

bfal_ipsl_200m_ts <- ts(bfal_ipsl_200m_ts_raw$Prob, start = c(1980,1), end = c(2100, 365), frequency = 365)
plot.ts(bfal_ipsl_200m_ts)

bfal_ipsl_200m_ts_decomp <- decompose(bfal_ipsl_200m_ts)
# print(bfal_ipsl_200m_ts_decomp)
plot(bfal_ipsl_200m_ts_decomp) # manually examine

#  bfal shelf-slope example ### !!!!!!!!!!!!!!
ipsl.dir = "bfal/dailyP_IPSL/"  # change to the folder where your daily projections are stored
setwd(ipsl.dir)
list.files()

bfal_ipsl_shelfslope_ts_raw <- read.csv("bfal43_MeanOccProb_TimeSeries_1980-2100_ShelfSlope.csv")
head(bfal_ipsl_shelfslope_ts_raw)
bfal_ipsl_shelfslope_ts_raw$Date <- as.Date(bfal_ipsl_shelfslope_ts_raw$Date)
bfal_ipsl_shelfslope_ts_raw <- bfal_ipsl_shelfslope_ts_raw[!(month(bfal_ipsl_shelfslope_ts_raw$Date) == 2 & day(bfal_ipsl_shelfslope_ts_raw$Date) == 29), ] # remove leap day obs
bfal_ipsl_shelfslope_ts_raw <- bfal_ipsl_shelfslope_ts_raw[,-1]
bfal_ipsl_shelfslope_ts_raw$Mod <- 'IPSL'
bfal_ipsl_shelfslope_ts_raw$Area <- 'ShelfSlope'
bfal_ipsl_shelfslope_ts_raw$AreaType <- 'CCE'

bfal_ipsl_shelfslope_ts <- ts(bfal_ipsl_shelfslope_ts_raw$Prob, start = c(1980,1), end = c(2100, 365), frequency = 365)
plot.ts(bfal_ipsl_shelfslope_ts)

bfal_ipsl_shelfslope_ts_decomp <- decompose(bfal_ipsl_shelfslope_ts)
# print(bfal_ipsl_shelfslope_ts_decomp)
plot(bfal_ipsl_shelfslope_ts_decomp)


# __________________________ ####
# Calculate time of emergence (TOE) ####
data.dir <- "~/SeabirdModels/Data" # * CHANGE TO YOUR DATA DIRECTORY
setwd(data.dir)

all_cce_ts_raw <- read.csv('TimeSeriesAveragesByArea_AllSpecies.csv')

head(all_cce_ts_raw)

calculate_toe <- function(data) {
  # Step 2.1: Calculate baseline (mean and sd) from the historical data (1985-2015)
  baseline_data <- data %>% dplyr::filter(Date >= as.Date('1985-01-01') & Date <= as.Date('2015-12-31'))
  baseline_mean <- mean(baseline_data$Prob, na.rm = TRUE)
  baseline_sd <- sd(baseline_data$Prob, na.rm = TRUE)
  
  # Step 2.2: Define threshold (mean + 1 SD)
  threshold <- baseline_mean - (1 * baseline_sd)
  
  # Step 2.3: Identify the first time when Prob exceeds the threshold
  emergence_index <- which(data$Trend < threshold)[1]
  
  if (!is.na(emergence_index)) {
    emergence_year <- year(data$Date[emergence_index])
  } else {
    emergence_year <- NA  # If no emergence occurs (e.g., if the series doesn't exceed the threshold)
  }
  
  return(data.frame(Threshold = threshold, baseline_mean = baseline_mean, baseline_sd = baseline_sd, EmergenceYear = emergence_year))
}

# Step 3: Apply the function to each Species * Area * Mod group
toes_allcce <- all_cce_ts_raw %>%
  group_by(Species, Mod, Area) %>%
  do(calculate_toe(.)) %>%
  ungroup()

toes_allcce$EmergenceYearDate <- as.Date(paste(toes_allcce$EmergenceYear, "01", "01", sep = "-"), format = "%Y-%m-%d")

# View the results
print(toes_allcce, n = 150)
head(toes_allcce)

# Plot TOE ####
all_cce_ts_raw_with_toe <- all_cce_ts_raw %>%
  left_join(toes_allcce, by = c("Species", "Area", "Mod"))

unique(all_cce_ts_raw_with_toe$Mod)
head(all_cce_ts_raw_with_toe)

all_cce_ts_raw_with_toe$Date <- as.Date(all_cce_ts_raw_with_toe$Date)
all_cce_ts_raw_with_toe$EmergenceYearDate <- as.Date(all_cce_ts_raw_with_toe$EmergenceYearDate)


# Example plot:
ggplot() + 
  geom_line(data = subset(all_cce_ts_raw_with_toe, Area == 'CHNMS'), # using CHNMS as an example to plot
            aes(x = Date, y = Trend, col = Mod), alpha = 0.8) +
  geom_vline(data = subset(all_cce_ts_raw_with_toe, Area == 'CHNMS'), 
             aes(col = Mod, xintercept = EmergenceYearDate), size = 0.75, alpha = 0.75) +
  facet_wrap(~ Species, ncol = 5) +
  theme_Publication() + 
  scale_colour_Publication() +
  labs(y = 'HSI Trend') +
  theme(legend.position = "top", 
        legend.spacing = unit(0.1,'cm'),
        legend.key.size = unit(0.5,'cm'),
        axis.title.x = element_blank())


