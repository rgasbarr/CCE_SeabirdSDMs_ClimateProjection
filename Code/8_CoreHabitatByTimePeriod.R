# HEADER ####
# Author: Ryan Gasbarro; rygasbar@ucsc.edu
# Script Purpose: Calculate changes in projected habitat suitability for seabird species by time period (e.g. historical vs. future)
# CCE seabird species using habitat suitability model outputs
# See Gasbarro et al. (2025) PLoS Clim. for information
# PACKAGES & FUNCITONS ####
library(glue)
library(geosphere)
library(tidyverse)
library(wesanderson)

theme_Publication <- function(base_size=12, base_family="helvetica") {
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
            legend.key.size= unit(0.4, "cm"),
            legend.spacing = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            legend.box.background = element_rect(size = 1),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_blank(),
            strip.text = element_text(size = 10)
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

# LOAD DATA & CALCULATE ####
# _____________
data.dir <- "~/SeabirdModels/Data" # ***** CHANGE TO YOUR DATA DIRECTORY *****
setwd(data.dir)

# NOTE: All time-series must be created in step 4 and combined [e.g. using rbind()] before completing this step **********************
# Alternatively, all combined time-series are provided as 'TimeSeriesAveragesByArea_AllSpecies.csv

all_cce_ts_raw <-  read.csv('TimeSeriesAveragesByArea_AllSpecies.csv')[,-1] # delete first 2 column: non-meaningful data

# all_cce_ts_caau <- subset(all_cce_ts_raw, Species == 'CAAU' & Area == 'ShelfSlope') # this Area is the 273 km offshore area covered by both ship- and air-based surveys
# all_cce_ts_rhau <- subset(all_cce_ts_raw, Species == 'RHAU' & Area == 'ShelfSlope')
# all_cce_ts_comu <- subset(all_cce_ts_raw, Species == 'COMU' & Area == 'ShelfSlope')
# all_cce_ts_bfal <- subset(all_cce_ts_raw, Species == 'BFAL' & Area == 'ShelfSlope')
# all_cce_ts_sosh <- subset(all_cce_ts_raw, Species == 'SOSH' & Area == 'ShelfSlope')

all_cce_ts_raw$doy <- yday(all_cce_ts_raw$Date) # calculate Julian day

period_cat <- function(years) {
  ifelse(years >= 1985 & years <= 2015, "1985-2015",
         ifelse(years >= 2035 & years <= 2065, "2035-2065",
                ifelse(years >= 2070 & years <= 2100, "2070-2100",
                       NA)))
}


# SUMMARIZE BY PERIOD
summary_cce_ts <- all_cce_ts_raw %>% dplyr::group_by(Species, period, doy) %>% 
  dplyr::summarize(avg_prob = mean(Prob, na.rm = TRUE), # avg pr
                   sd_prob = sd(Prob, na.rm = TRUE),
                   avg_trend = mean(Trend, na.rm = TRUE),
                   sd_trend = sd(Trend, na.rm = TRUE),
                   avg_rand = mean(Rand, na.rm = TRUE),
                   sd_rand = sd(Rand, na.rm = TRUE),
                   avg_seas = mean(Seas, na.rm = TRUE),
                   sd_seas = sd(Seas, na.rm = TRUE))

head(summary_cce_ts)


pal12 <- wes_palette("Zissou1", 12, type = "continuous")

summary_cce_ts_period <- summary_cce_ts %>%
  filter(!is.na(period))

summary_cce_ts_period <- summary_cce_ts_period %>%
  mutate(period2 = case_when(
    period == "1985-2015" ~ 2000,
    period == "2035-2065" ~ 2050,
    period == "2070-2100" ~ 2085,
    TRUE ~ NA_real_  # In case there are other values
  ))


# calculate prop. of days with average HSI over core habitat threshold for each area for each time period
species_labels <- c(
  "BFAL" = "Black-footed albatross",
  "CAAU" = "Cassin's auklet",
  "COMU" = "Common murre",
  "RHAU" = "Rhinoceros auklet",
  "SOSH" = "Sooty shearwater"
)

# use historical averages over the study area (273 km offshore distance) to calculate 'core habitat'
# i.e. >75th percentile of historical HSI
shelfslope_filt <- all_cce_ts_raw %>% 
  filter(period == '1985-2015', Area == 'ShelfSlope')

per75_valu <- shelfslope_filt %>%
  group_by(Species) %>%
  dplyr::summarize(per75 = quantile(Prob, probs = 0.75, na.rm = TRUE))

all_cce_ts_raw <- all_cce_ts_raw %>%
  left_join(per75_valu, by = "Species")

# manually inspect
# unique(all_cce_ts_raw$Area) 
# summary(as.factor(all_cce_ts_raw$Area))

all_cce_ts_raw$Area[all_cce_ts_raw$Area == 'ShelfSlope'] <- 'Study.Area' # optional rename

# add province variable for plotting
all_cce_ts_raw <- 
  all_cce_ts_raw %>%
  mutate(province = case_when(
    Area == "200m" ~ "NA",
    Area == "Study.Area" ~ "NA", # optional rename
    Area == "MBNMS" ~ "Central",
    Area == "GFNMS" ~ "Central",
    Area == "CBNMS" ~ "Central",
    Area == "CHNMS" ~ "Central",
    Area == "MorroBay" ~ "Central",
    Area == "Humboldt" ~ "North",
    Area == "Brookings" ~ "North",
    Area == "CoosBay" ~ "North",
    TRUE ~ "NA"  # Default case
  ))

pal4 <- wes_palette("Moonrise2")

all_cce_ts_raw$province <- factor(all_cce_ts_raw$province, levels = c('South', 'Central', 'North', 'NA'))

all_cce_ts_raw$Area <- factor(all_cce_ts_raw$Area, levels = c('200m', 'Study.Area',
                                                              'MBNMS', 'GFNMS', 'CBNMS', 'CHNMS',
                                                              'MorroBay', 'Humboldt', 'Brookings', 'CoosBay'))



# summary(all_cce_ts_raw$Area) # manually inspect

all_cce_ts_raw2 <- all_cce_ts_raw %>%
  filter(!is.na(period))

# create df that contains proportion of days over core habitat threshold
all_cce_ts_dat <- all_cce_ts_raw2 %>% 
  group_by(Species,Area, period) %>%
  dplyr::summarize(percentage_above_per75 = mean(Prob > per75, na.rm = TRUE) * 100)

# add province variable to this df
all_cce_ts_dat <- 
  all_cce_ts_dat %>%
  mutate(province = case_when(
    Area == "200m" ~ "NA",
    Area == "Study.Area" ~ "NA", # optional rename
    Area == "MBNMS" ~ "Central",
    Area == "GFNMS" ~ "Central",
    Area == "CBNMS" ~ "Central",
    Area == "CHNMS" ~ "Central",
    Area == "MorroBay" ~ "Central",
    Area == "Humboldt" ~ "North",
    Area == "Brookings" ~ "North",
    Area == "CoosBay" ~ "North",
    TRUE ~ "NA"  # Default case
  ))

# PLOT ####
pals6 <- wes_palette('AsteroidCity1')

all_cce_ts_dat$Species <- factor(all_cce_ts_dat$Species, levels = c('BFAL', 'SOSH', 'COMU', 'CAAU', 'RHAU'))

ggplot(all_cce_ts_dat, aes(y = Area, x = percentage_above_per75, fill = period)) +
  geom_bar(stat = "identity", position = "identity", size = 2) +
  scale_fill_manual(values = pals6) +  # Use fill for color
  facet_wrap(~Species, ncol = 5, labeller = labeller(Species = species_labels)) +
  labs(
    y = "Area", 
    x = "% Days", 
    fill = "Period",
    alpha = 'Hist. Suitability') +
  theme_Publication() + 
  geom_hline(yintercept = c(4.5, 8.5), linetype = "dashed", color = "black") +  # Add dashed lines
  scale_y_discrete(limits = rev(levels(all_cce_ts_dat$Area))) +  
  theme(
    plot.margin = unit(c(0.25, 0.5, 0.15, 0.25), "cm"),
    legend.background = element_blank(),
    legend.box.background = element_rect(size = 0.33), 
    strip.text = element_text(size = 14))

# Make ensemble mean tileplots for each area
esm_mean_ts_all <- all_cce_ts_raw2 %>%
  group_by(Species, Area, doy, period) %>%
  dplyr::summarise(esmMean = mean(Prob, na.rm = TRUE)) %>%
  ungroup()  # Ungroup after summarising

esm_mean_ts_all <- esm_mean_ts_all %>%
  group_by(Species, Area, doy, period) %>%
  mutate(
    qrt75 = quantile(esmMean, probs = 0.75), # Calculate 75th perc
    above_threshold = esmMean > qrt75  # Check if esmMean is above
  )


zissou_colors <- wes_palette(name = "Zissou1", n = 100, type = 'continuous')

esm_mean_ts_all$Species <- factor(esm_mean_ts_all$Species, 
                              levels = c('BFAL', 'SOSH', 'COMU', 'CAAU', 'RHAU'),
                              labels = c('BFAL', 'SOSH', 'COMU', 'CAAU', 'RHAU'))

levels(esm_mean_ts_all$Species)

esm_mean_ts_all %>%
  ggplot(aes(y = period, x = doy)) +
  geom_raster(aes(fill = esmMean)) +
  geom_point(aes(col = above_threshold), size = 2) +  # Set size for points
  facet_grid(rows = vars(Area), cols = vars(Species)) + 
  scale_fill_gradientn("Suitability", colours = zissou_colors, na.value = "black") +
  geom_vline(aes(xintercept = 73), col = 'black',  lty = 2) +
  geom_vline(aes(xintercept = 227), col = 'black',  lty = 2) +
  geom_vline(aes(xintercept = 319), col = 'black',  lty = 2) +
  scale_y_discrete(expand = c(0, 0)) +  # Use scale_y_discrete()
  scale_x_continuous(expand = c(0, 0)) +
  xlab("Day") +
  ylab("Period") +
  theme_Publication() +
  theme(legend.position = "right", 
        legend.direction = 'vertical', 
        strip.placement = 'inside',
        strip.background = element_blank(),  
        panel.border = element_rect(color = "black", size = 1),
        panel.spacing = unit(0.75, "mm")) +
  guides(col = 'none') +
  scale_color_manual(values = c("TRUE" = "white", "FALSE" = 'transparent')) 


