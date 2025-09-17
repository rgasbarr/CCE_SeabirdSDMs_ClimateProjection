# HEADER ####
# Author: Ryan Gasbarro; rygasbar@ucsc.edu
# Purpose: Plot presence/absence data for five seabird species in the California Current
# See Gasbarro et al. (2025) PLoS Climate for more information

# LOAD PACKAGES ####
library(tidyverse)
library(gridExtra)
library(sf)
library(ggsci)
library(ggthemes)
library(adehabitatHR)
library(rgeos)
library(maptools)
library(sp)
library(rnaturalearth)
library(rnaturalearthdata)
library(wesanderson)
library(ggpubfigs)
library(ggpubr)
library(lubridate)


# LOAD & CLEAN DATA ####
rm(list = ls())
set.seed(42)
data.dir <- "~/SeabirdModels/Data" # * CHANGE TO YOUR DATA DIRECTORY
setwd(data.dir)
dat <- read.csv("SeabirdsExampleData.csv") # df with seabird obs, metadata (e.g. date, lat, long, observing platform), matched ROMS (e.g. temperature) & static (e.g. bathymetry, distance-to-shore) data
dat <- dat[,-1] # remove unneccessary col

# create a factor for biogeographic province based on latitudinal breaks at Pt. Conception & Mendocino
dat$province <- factor("")
dat$province <- ifelse(dat$lat > 40.4401, "NORTH", dat$province)
dat$province <- ifelse(dat$lat < 40.4401 & dat$lat > 34.4486, "CENTRAL", dat$province)
dat$province <- ifelse(dat$lat < 34.4486, "SOUTH", dat$province)
dat$province <- as.factor(dat$province)

# Create a function to assign oceanographic season
get_season <- function(date) {
  month <- month(date)
  day <- day(date)
  ifelse((month == 3 & day >= 15) | (month > 3 & month < 8) | (month == 8 & day <= 14), "Upwelling",
         ifelse((month == 8 & day >= 15) | (month > 8 & month < 11) | (month == 11 & day <= 14), "Oceanic",
                "Davidson"))
}

dat <- dat %>% mutate(season = factor(get_season(date)))

# 134, -115.5, 30, 48  (xmin, xmax, ymin, ymax)
dat <- subset(dat, lat > 30 & lat <48 & lon > -134 & lon < -115.5) # remove points not on ROMS grid
dat$jday <- lubridate::yday(dat$date)

dat_boat <- subset(dat, platform == 'boat')
dat_aerial <- subset(dat, platform == 'aerial')

max(dat_aerial$distance_to_shore) # ~273000 is the max. distance the aerial surveys go offshore

# Select the columns with spp counts
count_cols1 <- dat[, 15:19]

# Convert count variables to presence-absence format
presence_absence_cols1 <- ifelse(count_cols1 > 0, 1, 0)

dat[, 15:19] <- presence_absence_cols1

dat_pa <- gather(dat, key = "Species", value = "presAbs", 15:19)
dat_pa_air <- subset(dat_pa, platform == 'aerial')
dat_pa_boat <- subset(dat_pa, platform == 'boat')

dat_pa_air$alpha <- ifelse(dat_pa_air$presAbs == 0, 0.0000001, 1) # create a col to set ggplot alpha so absences don't swamp the pot
dat_pa_boat$alpha <- ifelse(dat_pa_boat$presAbs == 0, 0.0000001, 1) # create a col to set ggplot alpha so absences don't swamp the pot

dat_pa$season <- ordered(dat_pa$season, levels = c("Davidson", "Upwelling", "Oceanic")) # order the season var

setwd(data.dir)
z <- as.data.frame(rasterToPoints(raster('z.grd'))) # load bathymetry raster data for plotting
z$z2 <- round(z$layer) # round for cleaner plot
names(z) <- c('lon', 'lat', 'z','z_rounded')
head(z) # inspect

coastline <- ne_countries(scale = "medium", returnclass = 'sf') # load coastline obj
coastline_bbox <- c(xmin=-134,xmax=-116,ymin=30,ymax=50) # set boundaries
coastline_sf <- coastline %>% filter(admin %in% c("Canada", "Mexico", "United States of America")) %>% st_crop(coastline_bbox) # limit to relevant countries

states <- ne_states(country = "united states of america", returnclass = "sf") # load states obj
states <- states %>%filter(name %in% c('Washington', 'Oregon', 'California', 'Nevada', 'Idaho')) %>%  # filter relevant states
  st_crop(coastline_bbox)

shape.dir <- '~/SeabirdModels/Shapefiles' # directory where shapefiles are stored
setwd(shape.dir)
sanc <- read_sf('CCE_NMS_combined_shapefile.shp') # contains all sanctuaries except Chumash
chnms <- read_sf('CHNMS_Boundary_10152024') #  Chumash

wind.dir <- '~/SeabirdModels/Shapefiles/WindAreas' # diretory where wind area shapefiles are stored
ca <- read_sf(paste(wind.dir, '/CA_Wind_Shapefiles/LeaseAreas_CA_Outlines_2022_09_08.shp', sep = '')) # two CA areas
or <- read_sf(paste(wind.dir, '/OR_Wind_Shapefiles/OREGON_Wind_Energy_Proposed_Lease_Areas_2024_BOEM.shp', sep = '')) # two OR areas

# PLOT ####

# theme for plotting
my.theme <- theme(
  panel.background = element_blank(),
  panel.border = element_rect(colour="black", fill = NA),
  panel.grid = element_blank(),
  axis.line = element_line("black"),
  text = element_text(size = 11),
  axis.text = element_text(size = 11, colour = "black"),
  axis.title = element_text(size = 12, colour = "black"),
  legend.title = element_text(size = 12),
  legend.text = element_text(size = 11),
  legend.key = element_blank(),
  legend.background = element_blank(),
  legend.position = 'top')

# Plot boat surveys
surv_boat <- ggplot() + 
  geom_tile(data = subset(z, z < 0 & dist < 273000), aes(x = lon, y = lat, fill = z_rounded)) +
  geom_point(data = subset(dat_pa, Species == 'sosh' & platform == 'boat' & distance_to_shore < 273000), aes(x = lon, y = lat), fill = 'grey80', col = 'grey80', size = 0.4, alpha = 0.1) + 
  geom_sf(data = coastline_sf, fill = 'transparent', col = 'black') + 
  geom_sf(data = states, fill = 'white', color = 'black') +
  geom_sf(data = sanc, col = 'black', fill = 'transparent', lwd = 0.52) +
  geom_sf(data = chnms, col = 'black', fill = 'transparent', lwd = 0.52) +
  geom_sf(data = ca, fill = 'white', color = 'black', lwd = 0.52, alpha = 0.8) +
  geom_sf(data = or, fill = 'white', color = 'black', lwd = 0.52, alpha = 0.8) +
  scale_x_continuous(limits = c(-128, -117), expand = c(0,0), breaks = c(-127, -124, -121, -118)) +
  scale_y_continuous(limits = c(32, 48.6), expand = c(0,0)) +
  labs(x = NULL, y = NULL, fill = 'Depth [m]', shape = 'Platform') +
  guides(alpha='none',
         # color = guide_legend(override.aes = list(size = 2, alpha = 1)),
         shape = guide_legend(override.aes = list(size = 2, alpha = 1, col = 'black'))) + 
  my.theme + 
  theme(legend.position = 'inside',
        legend.position.inside = c(0.8, 0.85), 
        legend.background = element_rect(fill = 'white', colour = 1),
        legend.direction = 'vertical',
        strip.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.text.y = element_blank(),
        axis.line = element_blank(),
        panel.spacing = unit(0, "cm"),
        plot.margin = margin(3, 0, 0, 0),
        strip.text.x = element_blank())

# Plot air surveys
surv_air <- ggplot() + 
  geom_tile(data = subset(z, z < 0 & dist < 273000), aes(x = lon, y = lat, fill = z_rounded)) +
  geom_point(data = subset(dat_pa, Species == 'sosh' & platform == 'aerial'), aes(x = lon, y = lat), fill = 'grey80', col = 'grey80', size = 0.4, alpha = 0.1) + 
  geom_sf(data = coastline_sf, fill = 'transparent', col = 'black') + 
  geom_sf(data = states, fill = 'white', color = 'black') +
  geom_sf(data = sanc, col = 'black', fill = 'transparent', lwd = 0.52) +
  geom_sf(data = chnms, col = 'black', fill = 'transparent', lwd = 0.52) +
  geom_sf(data = ca, fill = 'white', color = 'black', lwd = 0.52, alpha = 0.8) +
  geom_sf(data = or, fill = 'white', color = 'black', lwd = 0.52, alpha = 0.8) +
  scale_x_continuous(limits = c(-128, -117), expand = c(0,0), breaks = c(-127, -124, -121, -118)) +
  scale_y_continuous(limits = c(32, 48.6), expand = c(0,0)) +
  labs(x = NULL, y = NULL, fill = 'Depth [m]', shape = 'Platform') +
  guides(alpha='none', fill = 'none',
         # color = guide_legend(override.aes = list(size = 2, alpha = 1)),
         shape = guide_legend(override.aes = list(size = 2, alpha = 1, col = 'black'))) + 
  my.theme + 
  theme(legend.position = 'inside',
        legend.position.inside = c(0.8, 0.82), 
        legend.background = element_rect(fill = 'white', colour = 1),
        legend.direction = 'vertical',
        strip.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.text.y = element_blank(),
        panel.spacing = unit(0, "cm"),
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = margin(3, 0, 0, 0),
        strip.text.x = element_blank())


# Plots of presences in each ocean season by species ## !!!
bfal_pres <- ggplot() + 
  geom_point(data = subset(dat_pa, Species == 'bfal' & presAbs == 1 & distance_to_shore < 273000), aes(x = lon, y = lat, col = season),  size = 1,  alpha = 0.1) + 
  geom_sf(data = coastline_sf, fill = 'transparent', color = 'black') + 
  geom_sf(data = states, fill = 'white', color = 'black') +
  # annotate('text', label = italic('Black-footed\nalbatross'), x = -119.5, y= 47, size = 6) +
  geom_sf(data = sanc, fill = NA, color = 'black', lwd = 0.52) +
  geom_sf(data = chnms, fill = NA, color = 'black', lwd = 0.52) +
  geom_sf(data = ca, fill = 'white', color = 'black', lwd = 0.52, alpha = 0.8) +
  geom_sf(data = or, fill = 'white', color = 'black', lwd = 0.52, alpha = 0.8) +
  scale_x_continuous(limits = c(-128, -117), expand = c(0,0), breaks = c(-127, -124, -121, -118)) +
  scale_y_continuous(limits = c(32, 48.6), expand = c(0,0)) +
  scale_color_manual(values = friendly_pal("contrast_three")) +  
  labs(x = NULL, y = NULL, shape = 'Platform', col = 'Season') +
  guides(alpha='none', shape = 'none',
         color = guide_legend(override.aes = list(size = 2, alpha = 1))) + 
  my.theme + 
  theme(legend.position = 'inside',
        axis.text.y = element_blank(),
        legend.position.inside = c(0.7, 0.75), 
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        plot.margin = margin(3, 0, 0, 0),
        strip.background = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line = element_blank(),
        legend.background = element_rect(fill = 'white', colour = 1),
        panel.spacing = unit(0, "cm"),
        legend.direction = 'vertical')

caau_pres <- ggplot() + 
  geom_point(data = subset(dat_pa, Species == 'caau' & presAbs == 1 & distance_to_shore < 273000), aes(x = lon, y = lat, col = season),  size = 1, alpha = 0.1) + 
  geom_sf(data = coastline_sf, fill = 'transparent', color = 'black') + 
  geom_sf(data = states, fill = 'white', color = 'black') +
  # annotate('text', label = 'Cassin\'s\nauklet', x = -119.5, y= 47, size = 6) +
  geom_sf(data = sanc, fill = NA, color = 'black', lwd = 0.52) +
  geom_sf(data = chnms, fill = NA, color = 'black', lwd = 0.52) +
  geom_sf(data = ca, fill = 'white', color = 'black', lwd = 0.52, alpha = 0.8) +
  geom_sf(data = or, fill = 'white', color = 'black', lwd = 0.52, alpha = 0.8) +
  scale_x_continuous(limits = c(-128, -117), expand = c(0,0), breaks = c(-127, -124, -121, -118)) +
  scale_y_continuous(limits = c(32, 48.6), expand = c(0,0)) +
  scale_color_manual(values = friendly_pal("contrast_three")) +  
  labs(x = NULL, y = NULL, shape = 'Platform', col = 'Season') +
  guides(alpha='none', shape = 'none', col = 'none') + 
  my.theme + theme(axis.text.y = element_blank(), 
                   axis.ticks.y = element_blank(),
                   axis.line = element_blank(),
                   strip.background = element_blank(),
                   panel.border = element_rect(color = "black", fill = NA, size = 1), 
                   panel.spacing = unit(0, "cm"),
                   plot.margin = margin(3, 0, 0, 0))

comu_pres <- ggplot() + 
  geom_point(data = subset(dat_pa, Species == 'comu' & presAbs == 1 & distance_to_shore < 273000), aes(x = lon, y = lat, col = season),  size = 1, alpha = 0.1) + 
  geom_sf(data = coastline_sf, fill = 'transparent', color = 'black') + 
  geom_sf(data = states, fill = 'white', color = 'black') +
  # annotate('text', label = 'Common\nmurre', x = -119.5, y= 47, size = 6) +
  geom_sf(data = sanc, fill = NA, color = 'black', lwd = 0.52) +
  geom_sf(data = chnms, fill = NA, color = 'black', lwd = 0.52) +
  geom_sf(data = ca, fill = 'white', color = 'black', lwd = 0.52, alpha = 0.8) +
  geom_sf(data = or, fill = 'white', color = 'black', lwd = 0.52, alpha = 0.8) +
  scale_x_continuous(limits = c(-128, -117), expand = c(0,0), breaks = c(-127, -124, -121, -118)) +
  scale_y_continuous(limits = c(32, 48.6), expand = c(0,0)) +
  scale_color_manual(values = friendly_pal("contrast_three")) +  
  labs(x = NULL, y = NULL, shape = 'Platform', col = 'Season') +
  guides(alpha='none', shape = 'none', col = 'none') + 
  my.theme + theme(axis.text.y = element_blank(), 
                   axis.ticks.y = element_blank(),
                   axis.line = element_blank(),
                   strip.background = element_blank(),
                   panel.border = element_rect(color = "black", fill = NA, size = 1),
                   panel.spacing = unit(0, "cm"),
                   plot.margin = margin(3, 0, 0, 0))

rhau_pres <- ggplot() + 
  geom_point(data = subset(dat_pa, Species == 'rhau' & presAbs == 1 & distance_to_shore < 273000), aes(x = lon, y = lat, col = season),  size = 1, alpha = 0.1) + 
  geom_sf(data = coastline_sf, fill = 'transparent', color = 'black') + 
  geom_sf(data = states, fill = 'white', color = 'black') +
  # annotate('text', label = 'Rhinoceros\nauklet', x = -119.5, y= 47, size = 6) +
  geom_sf(data = sanc, fill = NA, color = 'black', lwd = 0.52) +
  geom_sf(data = chnms, fill = NA, color = 'black', lwd = 0.52) +
  geom_sf(data = ca, fill = 'white', color = 'black', lwd = 0.52, alpha = 0.8) +
  geom_sf(data = or, fill = 'white', color = 'black', lwd = 0.52, alpha = 0.8) +
  scale_x_continuous(limits = c(-128, -117), expand = c(0,0), breaks = c(-127, -124, -121, -118)) +
  scale_y_continuous(limits = c(32, 48.6), expand = c(0,0)) +
  scale_color_manual(values = friendly_pal("contrast_three")) +  
  labs(x = NULL, y = NULL, shape = 'Platform', col = 'Season') +
  guides(alpha='none', shape = 'none', col = 'none') + 
  my.theme + 
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.line = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        panel.spacing = unit(0, "cm"),
        plot.margin = margin(3, 0, 0, 0))

sosh_pres <- ggplot() + 
  geom_point(data = subset(dat_pa, Species == 'sosh' & presAbs == 1 & distance_to_shore < 273000), aes(x = lon, y = lat, col = season),  size = 1, alpha = 0.1) + 
  geom_sf(data = coastline_sf, fill = 'transparent', color = 'black') + 
  geom_sf(data = states, fill = 'white', color = 'black') +
  # annotate('text', label = 'Sooty\nshearwater', x = -119.5, y= 47, size = 6) +
  geom_sf(data = sanc, fill = NA, color = 'black', lwd = 0.52) +
  geom_sf(data = chnms, fill = NA, color = 'black', lwd = 0.52) +
  geom_sf(data = ca, fill = 'white', color = 'black', lwd = 0.52, alpha = 0.8) +
  geom_sf(data = or, fill = 'white', color = 'black', lwd = 0.52, alpha = 0.8) +
  scale_x_continuous(limits = c(-128, -117), expand = c(0,0), breaks = c(-127, -124, -121, -118)) +
  scale_y_continuous(limits = c(32, 48.6), expand = c(0,0)) +
  scale_color_manual(values = friendly_pal("contrast_three")) +  
  labs(x = NULL, y = NULL, shape = 'Platform', col = 'Season') +
  guides(alpha='none', shape = 'none', col = 'none') + 
  my.theme + theme(axis.text.y = element_blank(), 
                   axis.ticks.y = element_blank(),
                   axis.line = element_blank(),
                   strip.background = element_blank(),
                   panel.border = element_rect(color = "black", fill = NA, size = 1),
                   panel.spacing = unit(0, "cm"),
                   plot.margin = margin(3, 1, 0, 0))

grid.arrange(
  surv_boat, surv_air,
             bfal_pres, sosh_pres, comu_pres, caau_pres, rhau_pres, ncol = 5, padding = unit(0.25, "mm"))



############################# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
############################# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
############################# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
############################# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
dat_pa_air$Species <- factor(dat_pa_air$Species, 
                     levels = c('bfal', 'sosh', 'comu', 'caau', 'rhau'),
                     labels = c('BFAL', 'SOSH', 'COMU', 'CAAU', 'RHAU'))


dat_pa_boat$Species <- factor(dat_pa_boat$Species, 
                             levels = c('bfal', 'sosh', 'comu', 'caau', 'rhau'),
                             labels = c('BFAL', 'SOSH', 'COMU', 'CAAU', 'RHAU'))

# Check the changes
levels(dat_pa_air$Species)
levels(dat_pa_boat$Species)

# plot all species aerial surveys
all_sp_air_p <- ggplot() + 
  geom_point(data = subset(dat_pa_air, presAbs == 0), aes(x = lon, y = lat), color = 'lightgray', alpha = 0.1,  size = 0.3) + 
  geom_point(data = subset(dat_pa_air, presAbs == 1), aes(x = lon, y = lat, color = season),  size = 0.3, alpha = 0.3) + 
  facet_wrap(~Species, ncol = 5) + 
  geom_sf(data = coastline_sf, fill = 'transparent', color = 'black') + 
  geom_sf(data = states, fill = 'white', color = 'black') +
  # annotate('text', label = 'Sooty\nshearwater', x = -119.5, y= 47, size = 6) +
  geom_sf(data = sanc, fill = NA, color = 'black', lwd = 0.52) +
  geom_sf(data = chnms, fill = NA, color = 'black', lwd = 0.52) +
  geom_sf(data = ca, fill = 'white', color = 'black', lwd = 0.52, alpha = 0.8) +
  geom_sf(data = or, fill = 'white', color = 'black', lwd = 0.52, alpha = 0.8) +
  xlim(-126, -117.5) + ylim(33, 48) +  
  scale_color_manual(values = friendly_pal("contrast_three")) +  
  labs(x = NULL, y = NULL) +
  guides(alpha='none',
         color = guide_legend(override.aes = list(size = 3, alpha = 1))) + 
  my.theme + theme(legend.position = 'right', legend.direction = 'vertical', axis.text = element_text(size = 9))

# plot all spp ship-based surveys
all_sp_boat_p <- ggplot() + 
  geom_point(data = subset(dat_pa_boat, presAbs == 0), aes(x = lon, y = lat), color = 'lightgray', alpha = 0.1,  size = 0.3) + 
  geom_point(data = subset(dat_pa_boat, presAbs == 1), aes(x = lon, y = lat, color = season),  size = 0.3, alpha = 0.3) + 
  facet_wrap(~Species, ncol = 5) + 
  geom_sf(data = coastline_sf, fill = 'transparent', color = 'black') + 
  geom_sf(data = states, fill = 'white', color = 'black') +
  # annotate('text', label = 'Sooty\nshearwater', x = -119.5, y= 47, size = 6) +
  geom_sf(data = sanc, fill = NA, color = 'black', lwd = 0.52) +
  geom_sf(data = chnms, fill = NA, color = 'black', lwd = 0.52) +
  geom_sf(data = ca, fill = 'white', color = 'black', lwd = 0.52, alpha = 0.8) +
  geom_sf(data = or, fill = 'white', color = 'black', lwd = 0.52, alpha = 0.8) +
  xlim(-131, -117.5) + ylim(33, 48) +  
  scale_color_manual(values = friendly_pal("contrast_three")) +  
  labs(x = NULL, y = NULL) +
  guides(col = 'none', alpha='none') + my.theme + theme(axis.text = element_text(size = 9))

grid.arrange(all_sp_boat_p, all_sp_air_p)

