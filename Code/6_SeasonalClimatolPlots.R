# HEADER ####
# Author: Ryan Gasbarro; ryan.gasbarro@noaa.gov

# Script Purpose: Plot seasonal climatologies of occurrence probability from seabird SDMs

# *** NOTE *** Seasonal climatologies must be created before running this script
# See Script 4 ('4_CreateClimatologies&TimeSeries.R') M

# See Gasbarro et al. (2025) PLoS Clim. for information
# _______ #####
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
library(TTR)
library(gdata)
library(viridis)
library(viridisLite)
library(sf)
library(geosphere)
library(rmapshaper)
library(rnaturalearth)
library(rnaturalearthdata)


theme_Publication <- function(base_size=12, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation()
    + theme(plot.title = element_text(size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = 'black', size = 1),
            axis.title = element_text(size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(size = 12), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "top",
            legend.direction = "horizontal",
            legend.key.size= unit(0.4, "cm"),
            legend.spacing = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            legend.box.background = element_rect(size = 1),
            plot.margin=unit(c(2,0,0,0),"mm"),
            strip.background = element_blank(),
            strip.text = element_blank()
    ))
  
}

# Load coastline, sanctuaries, wind areas, etc. shapefiles ####
data.dir <- "~/SeabirdModels/Data" # * CHANGE TO YOUR DATA DIRECTORY
setwd(data.dir)

z <- as.data.frame(rasterToPoints(raster('z.grd')))
z_rast <- raster('z.grd')
z_rast <- round(z_rast)
z$z2 <- round(z$layer)
names(z) <- c('lon', 'lat', 'z','z_rounded')
head(z)

coastline <- ne_countries(scale = "medium", returnclass = 'sf')
coastline_bbox <- c(xmin=-134,xmax=-116,ymin=30,ymax=50)
coastline_sf <- coastline %>% 
  filter(admin %in% c("Canada", "Mexico", "United States of America")) %>% 
  st_crop(coastline_bbox)
dist <- raster('Distance_from_shore.grd')
dist_df <- as.data.frame(rasterToPoints(raster('Distance_from_shore.grd')))
states <- ne_states(country = "united states of america", returnclass = "sf")
states <- states %>%
  filter(name %in% c('Washington', 'Oregon', 'California', 'Nevada', 'Idaho')) %>% 
  st_crop(coastline_bbox)

# coastline_sf <- st_as_sf(coastline)
shape.dir <- '~/SeabirdModels/Shapefiles' # directory where National Marine Sanctuary (NMS) shapefiles are stored
sanc <- read_sf('Combined')
sanc$attribute1[4] <- 'CBNMS' 
chnms <- read_sf('ChumashHeritage')

wind.dir <- '~/SeabirdModels/Shapefiles/WindAreas' # directory where wind area shapefiles are stored
setwd(wind.dir)

ca <- read_sf('CA_Wind_Shapefiles/LeaseAreas_CA_Outlines_2022_09_08.shp')
or <- read_sf('OR_Wind_Shapefiles/OREGON_Wind_Energy_Proposed_Lease_Areas_2024_BOEM.shp')

# _______ #####
# Load Historical Seasonal Climatologies ####

seasonal_means_ipsl <- "~/SeabirdModels/dailyP_IPSL/SeasonalMeans"  # change to the folder where seasonal climatologies (see Script 4) are stored
seasonal_means_gfdl <- "~/SeabirdModels/dailyP_GFDL/SeasonalMeans"  # change to the folder where seasonal climatologies (see Script 4) are stored
seasonal_means_had <- "~/SeabirdModels/dailyP_HAD/SeasonalMeans"  # change to the folder where seasonal climatologies (see Script 4) are stored

# Load Rasters Davidson #
setwd(seasonal_means_gfdl)
bfal_gfdl_hist_davidson <- raster("Davidson_Historicalmean.tif")
caau_gfdl_hist_davidson <- raster("Davidson_Historicalmean.tif")
comu_gfdl_hist_davidson <- raster("Davidson_Historicalmean.tif")
rhau_gfdl_hist_davidson <- raster("Davidson_Historicalmean.tif")
sosh_gfdl_hist_davidson <- raster("Davidson_Historicalmean.tif")

setwd(seasonal_means_ipsl)
bfal_ipsl_hist_davidson <- raster("Davidson_Historicalmean.tif")
caau_ipsl_hist_davidson <- raster("Davidson_Historicalmean.tif")
comu_ipsl_hist_davidson <- raster("Davidson_Historicalmean.tif")
rhau_ipsl_hist_davidson <- raster("Davidson_Historicalmean.tif")
sosh_ipsl_hist_davidson <- raster("Davidson_Historicalmean.tif")

setwd(seasonal_means_had)
bfal_had_hist_davidson <- raster("Davidson_Historicalmean.tif")
caau_had_hist_davidson <- raster("Davidson_Historicalmean.tif")
comu_had_hist_davidson <- raster("Davidson_Historicalmean.tif")
rhau_had_hist_davidson <- raster("Davidson_Historicalmean.tif")
sosh_had_hist_davidson <- raster("Davidson_Historicalmean.tif")

# Load Rasters Upwelling #
setwd(seasonal_means_gfdl)
bfal_gfdl_hist_upwelling <- raster("Upwelling_Historicalmean.tif")
caau_gfdl_hist_upwelling <- raster("Upwelling_Historicalmean.tif")
comu_gfdl_hist_upwelling <- raster("Upwelling_Historicalmean.tif")
rhau_gfdl_hist_upwelling <- raster("Upwelling_Historicalmean.tif")
sosh_gfdl_hist_upwelling <- raster("Upwelling_Historicalmean.tif")

setwd(seasonal_means_ipsl)
bfal_ipsl_hist_upwelling <- raster("Upwelling_Historicalmean.tif")
caau_ipsl_hist_upwelling <- raster("Upwelling_Historicalmean.tif")
comu_ipsl_hist_upwelling <- raster("Upwelling_Historicalmean.tif")
rhau_ipsl_hist_upwelling <- raster("Upwelling_Historicalmean.tif")
sosh_ipsl_hist_upwelling <- raster("Upwelling_Historicalmean.tif")

setwd(seasonal_means_had)
bfal_had_hist_upwelling <- raster("Upwelling_Historicalmean.tif")
caau_had_hist_upwelling <- raster("Upwelling_Historicalmean.tif")
comu_had_hist_upwelling <- raster("Upwelling_Historicalmean.tif")
rhau_had_hist_upwelling <- raster("Upwelling_Historicalmean.tif")
sosh_had_hist_upwelling <- raster("Upwelling_Historicalmean.tif")

# Load Rasters Oceanic #
setwd(seasonal_means_gfdl)
bfal_gfdl_hist_oceanic <- raster("Oceanic_Historicalmean.tif")
caau_gfdl_hist_oceanic <- raster("Oceanic_Historicalmean.tif")
comu_gfdl_hist_oceanic <- raster("Oceanic_Historicalmean.tif")
rhau_gfdl_hist_oceanic <- raster("Oceanic_Historicalmean.tif")
sosh_gfdl_hist_oceanic <- raster("Oceanic_Historicalmean.tif")

setwd(seasonal_means_ipsl)
bfal_ipsl_hist_oceanic <- raster("Oceanic_Historicalmean.tif")
caau_ipsl_hist_oceanic <- raster("Oceanic_Historicalmean.tif")
comu_ipsl_hist_oceanic <- raster("Oceanic_Historicalmean.tif")
rhau_ipsl_hist_oceanic <- raster("Oceanic_Historicalmean.tif")
sosh_ipsl_hist_oceanic <- raster("Oceanic_Historicalmean.tif")

setwd(seasonal_means_had)
bfal_had_hist_oceanic <- raster("Oceanic_Historicalmean.tif")
caau_had_hist_oceanic <- raster("Oceanic_Historicalmean.tif")
comu_had_hist_oceanic <- raster("Oceanic_Historicalmean.tif")
rhau_had_hist_oceanic <- raster("Oceanic_Historicalmean.tif")
sosh_had_hist_oceanic <- raster("Oceanic_Historicalmean.tif")

# Plot GFDL ####
gfdl_seasonol_hist_stack <-raster::stack(bfal_gfdl_hist_davidson, bfal_gfdl_hist_upwelling, bfal_gfdl_hist_oceanic,
                                         caau_gfdl_hist_davidson, caau_gfdl_hist_upwelling, caau_gfdl_hist_oceanic,
                                         comu_gfdl_hist_davidson, comu_gfdl_hist_upwelling, comu_gfdl_hist_oceanic,
                                         rhau_gfdl_hist_davidson, rhau_gfdl_hist_upwelling, rhau_gfdl_hist_oceanic,
                                         sosh_gfdl_hist_davidson, sosh_gfdl_hist_upwelling, sosh_gfdl_hist_oceanic)

names(gfdl_seasonol_hist_stack) <- c('bfal_gfdl_hist_davidson', 'bfal_gfdl_hist_upwelling', 'bfal_gfdl_hist_oceanic',
                                     'caau_gfdl_hist_davidson', 'caau_gfdl_hist_upwelling', 'caau_gfdl_hist_oceanic',
                                     'comu_gfdl_hist_davidson', 'comu_gfdl_hist_upwelling', 'comu_gfdl_hist_oceanic',
                                     'rhau_gfdl_hist_davidson', 'rhau_gfdl_hist_upwelling', 'rhau_gfdl_hist_oceanic',
                                     'sosh_gfdl_hist_davidson', 'sosh_gfdl_hist_upwelling', 'sosh_gfdl_hist_oceanic')

gfdl_seasonol_hist_df <- as.data.frame(rasterToPoints(gfdl_seasonol_hist_stack, na.rm=FALSE))


gfdl_seasonol_hist_df <- gather(gfdl_seasonol_hist_df, model, prob, 3:17)
gfdl_seasonol_hist_df <- separate(data = gfdl_seasonol_hist_df, col = model, into = c("species", "season"), sep = "_gfdl_hist_")
gfdl_seasonol_hist_df$season <- factor(gfdl_seasonol_hist_df$season, levels = c("davidson", "upwelling","oceanic"))

# ggplot(gfdl_seasonol_hist_df) + geom_tile(aes(x=x,y=y, fill = prob)) + facet_wrap(~season+species, ncol = 5) +
#   geom_sf(data = coastline_sf, fill = 'transparent', col = 'black') +
#   geom_sf(data = states, fill = 'white', color = 'black') +
#   geom_sf(data = sanc, col = 'black', fill = 'transparent') +
#   geom_sf(data = ca, fill = 'transparent', color = 'black') +
#   geom_sf(data = or, fill = 'transparent', color = 'black') +
#   scale_x_continuous(limits = c(-128, -117), expand = c(0,0), breaks = c(-127, -124, -121, -118)) +
#   scale_y_continuous(limits = c(32, 47.25), expand = c(0,0)) +
#   scale_fill_viridis_c() +
#   guides(fill = 'none') +
#   # xlim(-126, -117.5) + ylim(30, 48) +  
#   labs(x = NULL, y = NULL) + theme_Publication()

# Plot IPSL ####
ipsl_seasonol_hist_stack <-raster::stack(bfal_ipsl_hist_davidson, bfal_ipsl_hist_upwelling, bfal_ipsl_hist_oceanic,
                                         caau_ipsl_hist_davidson, caau_ipsl_hist_upwelling, caau_ipsl_hist_oceanic,
                                         comu_ipsl_hist_davidson, comu_ipsl_hist_upwelling, comu_ipsl_hist_oceanic,
                                         rhau_ipsl_hist_davidson, rhau_ipsl_hist_upwelling, rhau_ipsl_hist_oceanic,
                                         sosh_ipsl_hist_davidson, sosh_ipsl_hist_upwelling, sosh_ipsl_hist_oceanic)

names(ipsl_seasonol_hist_stack) <- c('bfal_ipsl_hist_davidson', 'bfal_ipsl_hist_upwelling', 'bfal_ipsl_hist_oceanic',
                                     'caau_ipsl_hist_davidson', 'caau_ipsl_hist_upwelling', 'caau_ipsl_hist_oceanic',
                                     'comu_ipsl_hist_davidson', 'comu_ipsl_hist_upwelling', 'comu_ipsl_hist_oceanic',
                                     'rhau_ipsl_hist_davidson', 'rhau_ipsl_hist_upwelling', 'rhau_ipsl_hist_oceanic',
                                     'sosh_ipsl_hist_davidson', 'sosh_ipsl_hist_upwelling', 'sosh_ipsl_hist_oceanic')

ipsl_seasonol_hist_df <- as.data.frame(rasterToPoints(ipsl_seasonol_hist_stack, na.rm=FALSE))


ipsl_seasonol_hist_df <- gather(ipsl_seasonol_hist_df, model, prob, 3:17)
ipsl_seasonol_hist_df <- separate(data = ipsl_seasonol_hist_df, col = model, into = c("species", "season"), sep = "_ipsl_hist_")
ipsl_seasonol_hist_df$season <- factor(ipsl_seasonol_hist_df$season, levels = c("davidson", "upwelling","oceanic"))

# ggplot(ipsl_seasonol_hist_df) + geom_tile(aes(x=x,y=y, fill = prob)) + facet_wrap(~season+species, ncol = 5) +
#   geom_sf(data = coastline_sf, fill = 'black') +
#   geom_sf(data = states, fill = 'white', color = 'black') +
#   geom_sf(data = sanc, col = 'black', fill = 'transparent') +
#   geom_sf(data = ca, fill = 'transparent', color = 'black',  alpha = 0.75) +
#   geom_sf(data = or, fill = 'transparent', color = 'black',  alpha = 0.75) +
#   scale_x_continuous(limits = c(-128, -117), expand = c(0,0), breaks = c(-127, -124, -121, -118)) +
#   scale_y_continuous(limits = c(32, 47.25), expand = c(0,0)) +
#   scale_fill_viridis_c() +
#   guides(fill = 'none') +
#   # xlim(-126, -117.5) + ylim(30, 48) +  
#   labs(x = NULL, y = NULL) + theme_Publication()

# Plot HAD ####
had_seasonol_hist_stack <-raster::stack(bfal_had_hist_davidson, bfal_had_hist_upwelling, bfal_had_hist_oceanic,
                                         caau_had_hist_davidson, caau_had_hist_upwelling, caau_had_hist_oceanic,
                                         comu_had_hist_davidson, comu_had_hist_upwelling, comu_had_hist_oceanic,
                                         rhau_had_hist_davidson, rhau_had_hist_upwelling, rhau_had_hist_oceanic,
                                         sosh_had_hist_davidson, sosh_had_hist_upwelling, sosh_had_hist_oceanic)

names(had_seasonol_hist_stack) <- c('bfal_had_hist_davidson', 'bfal_had_hist_upwelling', 'bfal_had_hist_oceanic',
                                     'caau_had_hist_davidson', 'caau_had_hist_upwelling', 'caau_had_hist_oceanic',
                                     'comu_had_hist_davidson', 'comu_had_hist_upwelling', 'comu_had_hist_oceanic',
                                     'rhau_had_hist_davidson', 'rhau_had_hist_upwelling', 'rhau_had_hist_oceanic',
                                     'sosh_had_hist_davidson', 'sosh_had_hist_upwelling', 'sosh_had_hist_oceanic')

had_seasonol_hist_df <- as.data.frame(rasterToPoints(had_seasonol_hist_stack, na.rm=FALSE))


had_seasonol_hist_df <- gather(had_seasonol_hist_df, model, prob, 3:17)
had_seasonol_hist_df <- separate(data = had_seasonol_hist_df, col = model, into = c("species", "season"), sep = "_had_hist_")
had_seasonol_hist_df$season <- factor(had_seasonol_hist_df$season, levels = c("davidson", "upwelling","oceanic"))

# ggplot(had_seasonol_hist_df) + geom_tile(aes(x=x,y=y, fill = prob)) + facet_wrap(~season+species, ncol = 5) +
#   geom_sf(data = coastline_sf, fill = 'black') +
#   geom_sf(data = states, fill = 'white', color = 'black') +
#   geom_sf(data = sanc, col = 'black', fill = 'transparent') +
#   geom_sf(data = ca, fill = 'transparent', color = 'black') +
#   geom_sf(data = or, fill = 'transparent', color = 'black') +
#   scale_x_continuous(limits = c(-128, -117), expand = c(0,0), breaks = c(-127, -124, -121, -118)) +
#   scale_y_continuous(limits = c(32, 47.25), expand = c(0,0)) +
#   scale_fill_viridis_c() +
#   guides(fill = 'none') +
#   # xlim(-126, -117.5) + ylim(30, 48) +  
#   labs(x = NULL, y = NULL) + theme_Publication()

# Plot Ensemble ####
gfdl_seasonol_hist_df$mod <- 'gfdl'
ipsl_seasonol_hist_df$mod <- 'ipsl'
had_seasonol_hist_df$mod <- 'had'

ensem_hist_seasonclimat_df <- cbind(gfdl_seasonol_hist_df, ipsl_seasonol_hist_df[,c(5,6)], had_seasonol_hist_df[,c(5,6)])
ensem_hist_seasonclimat_df <- ensem_hist_seasonclimat_df[!duplicated(as.list(ensem_hist_seasonclimat_df))]

names(ensem_hist_seasonclimat_df)[5:10] <- c('prob_gfdl', 'mod', 'prob_ipsl', 'mod2', 'prob_had', 'mod3')
names(ensem_hist_seasonclimat_df)
ensem_hist_seasonclimat_df <- ensem_hist_seasonclimat_df[,-c(6, 8, 10)]

ensem_hist_seasonclimat_df$prob_em <- rowMeans(ensem_hist_seasonclimat_df[5:7], na.rm=TRUE) # calc ensemble mean
ensem_hist_seasonclimat_df$prob_esd <-apply(ensem_hist_seasonclimat_df[,5:7],1,sd)  # calc ensemble sd
head(ensem_hist_seasonclimat_df)
ensem_hist_seasonclimat_df$dist <- raster::extract(dist, ensem_hist_seasonclimat_df[,1:2])

#ensemble mean HSI plot 
ggplot(subset(ensem_hist_seasonclimat_df, dist < 273000)) + 
  geom_tile(aes(x=x,y=y, fill = prob_em)) + 
  facet_wrap(~season+species, ncol = 5) +
  geom_sf(data = coastline_sf, fill = 'black') +
  geom_sf(data = states, fill = 'white', color = 'black') +
  geom_sf(data = sanc, col = 'white', fill = 'transparent') +
  geom_sf(data = ca, fill = 'transparent', color = 'white') +
  geom_sf(data = or, fill = 'transparent', color = 'white') +
  scale_x_continuous(limits = c(-128, -117), expand = c(0,0), breaks = c(-127, -124, -121, -118)) +
  scale_y_continuous(limits = c(32, 47.25), expand = c(0,0)) +
  scale_fill_viridis_c() +
  # guides(fill = 'none') +
  labs(x = NULL, y = NULL, fill = 'HSI') + 
  theme_Publication() + theme(legend.position = c(0.95,0.95), 
                              legend.direction = 'vertical',
                              strip.background = element_blank())
#ensemble sd HSI plot 
ggplot(subset(ensem_hist_seasonclimat_df, dist < 273000)) + 
  geom_tile(aes(x=x,y=y, fill = prob_esd)) + 
  facet_wrap(~season+species, ncol = 5) +
  geom_sf(data = coastline_sf, fill = 'black') +
  geom_sf(data = states, fill = 'white', color = 'black') +
  geom_sf(data = sanc, col = 'white', fill = 'transparent') +
  geom_sf(data = ca, fill = 'transparent', color = 'white') +
  geom_sf(data = or, fill = 'transparent', color = 'white') +
  scale_x_continuous(limits = c(-128, -117), expand = c(0,0), breaks = c(-127, -124, -121, -118)) +
  scale_y_continuous(limits = c(32, 47.25), expand = c(0,0)) +
  scale_fill_viridis_c() +
  # guides(fill = 'none') +
  labs(x = NULL, y = NULL, fill = 'ESM SD') + 
  theme_Publication() + theme(legend.position = c(0.95,0.95), 
                              legend.direction = 'vertical', 
                              strip.background = element_blank(),
                              plot.margin = margin(1, 0, 0, 0))


# _______ #####
# Load Future Seasonal Climatologies ####
# Load Rasters Davidson #
seasonal_means_ipsl <- "~/SeabirdModels/dailyP_IPSL/SeasonalMeans"  # change to the folder where seasonal climatologies (see Script 4) are stored
seasonal_means_gfdl <- "~/SeabirdModels/dailyP_GFDL/SeasonalMeans"  # change to the folder where seasonal climatologies (see Script 4) are stored
seasonal_means_had <- "~/SeabirdModels/dailyP_HAD/SeasonalMeans"  # change to the folder where seasonal climatologies (see Script 4) are stored

setwd(seasonal_means_gfdl)
bfal_gfdl_fut_davidson <- raster("Davidson_Futuremean.tif")
caau_gfdl_fut_davidson <- raster("Davidson_Futuremean.tif")
comu_gfdl_fut_davidson <- raster("Davidson_Futuremean.tif")
rhau_gfdl_fut_davidson <- raster("Davidson_Futuremean.tif")
sosh_gfdl_fut_davidson <-raster("Davidson_Futuremean.tif")

setwd(seasonal_means_ipsl)
bfal_ipsl_fut_davidson <- raster("Davidson_Futuremean.tif")
caau_ipsl_fut_davidson <- raster("Davidson_Futuremean.tif")
comu_ipsl_fut_davidson <- raster("Davidson_Futuremean.tif")
rhau_ipsl_fut_davidson <- raster("Davidson_Futuremean.tif")
sosh_ipsl_fut_davidson <-raster("Davidson_Futuremean.tif")

setwd(seasonal_means_had)
bfal_had_fut_davidson <- raster("Davidson_Futuremean.tif")
caau_had_fut_davidson <- raster("Davidson_Futuremean.tif")
comu_had_fut_davidson <- raster("Davidson_Futuremean.tif")
rhau_had_fut_davidson <- raster("Davidson_Futuremean.tif")
sosh_had_fut_davidson <-raster("Davidson_Futuremean.tif")

# Load Rasters Upwelling #
setwd(seasonal_means_gfdl)
bfal_gfdl_fut_upwelling <- raster("Upwelling_Futuremean.tif")
caau_gfdl_fut_upwelling <- raster("Upwelling_Futuremean.tif")
comu_gfdl_fut_upwelling <- raster("Upwelling_Futuremean.tif")
rhau_gfdl_fut_upwelling <- raster("Upwelling_Futuremean.tif")
sosh_gfdl_fut_upwelling <-raster("Upwelling_Futuremean.tif")

setwd(seasonal_means_ipsl)
bfal_ipsl_fut_upwelling <- raster("Upwelling_Futuremean.tif")
caau_ipsl_fut_upwelling <- raster("Upwelling_Futuremean.tif")
comu_ipsl_fut_upwelling <- raster("Upwelling_Futuremean.tif")
rhau_ipsl_fut_upwelling <- raster("Upwelling_Futuremean.tif")
sosh_ipsl_fut_upwelling <-raster("Upwelling_Futuremean.tif")

setwd(seasonal_means_had)
bfal_had_fut_upwelling <- raster("Upwelling_Futuremean.tif")
caau_had_fut_upwelling <- raster("Upwelling_Futuremean.tif")
comu_had_fut_upwelling <- raster("Upwelling_Futuremean.tif")
rhau_had_fut_upwelling <- raster("Upwelling_Futuremean.tif")
sosh_had_fut_upwelling <-raster("Upwelling_Futuremean.tif")

# Load Rasters Oceanic #
setwd(seasonal_means_gfdl)
bfal_gfdl_fut_oceanic <- raster("Oceanic_Futuremean.tif")
caau_gfdl_fut_oceanic <- raster("Oceanic_Futuremean.tif")
comu_gfdl_fut_oceanic <- raster("Oceanic_Futuremean.tif")
rhau_gfdl_fut_oceanic <- raster("Oceanic_Futuremean.tif")
sosh_gfdl_fut_oceanic <-raster("Oceanic_Futuremean.tif")

setwd(seasonal_means_ipsl)
bfal_ipsl_fut_oceanic <- raster("Oceanic_Futuremean.tif")
caau_ipsl_fut_oceanic <- raster("Oceanic_Futuremean.tif")
comu_ipsl_fut_oceanic <- raster("Oceanic_Futuremean.tif")
rhau_ipsl_fut_oceanic <- raster("Oceanic_Futuremean.tif")
sosh_ipsl_fut_oceanic <-raster("Oceanic_Futuremean.tif")

setwd(seasonal_means_had)
bfal_had_fut_oceanic <- raster("Oceanic_Futuremean.tif")
caau_had_fut_oceanic <- raster("Oceanic_Futuremean.tif")
comu_had_fut_oceanic <- raster("Oceanic_Futuremean.tif")
rhau_had_fut_oceanic <- raster("Oceanic_Futuremean.tif")
sosh_had_fut_oceanic <-raster("Oceanic_Futuremean.tif")

# Plot GFDL ####
gfdl_seasonol_fut_stack <-raster::stack(bfal_gfdl_fut_davidson, bfal_gfdl_fut_upwelling, bfal_gfdl_fut_oceanic,
                                         caau_gfdl_fut_davidson, caau_gfdl_fut_upwelling, caau_gfdl_fut_oceanic,
                                         comu_gfdl_fut_davidson, comu_gfdl_fut_upwelling, comu_gfdl_fut_oceanic,
                                         rhau_gfdl_fut_davidson, rhau_gfdl_fut_upwelling, rhau_gfdl_fut_oceanic,
                                         sosh_gfdl_fut_davidson, sosh_gfdl_fut_upwelling, sosh_gfdl_fut_oceanic)

names(gfdl_seasonol_fut_stack) <- c('bfal_gfdl_fut_davidson', 'bfal_gfdl_fut_upwelling', 'bfal_gfdl_fut_oceanic',
                                     'caau_gfdl_fut_davidson', 'caau_gfdl_fut_upwelling', 'caau_gfdl_fut_oceanic',
                                     'comu_gfdl_fut_davidson', 'comu_gfdl_fut_upwelling', 'comu_gfdl_fut_oceanic',
                                     'rhau_gfdl_fut_davidson', 'rhau_gfdl_fut_upwelling', 'rhau_gfdl_fut_oceanic',
                                     'sosh_gfdl_fut_davidson', 'sosh_gfdl_fut_upwelling', 'sosh_gfdl_fut_oceanic')

gfdl_seasonol_fut_df <- as.data.frame(rasterToPoints(gfdl_seasonol_fut_stack, na.rm=FALSE))


gfdl_seasonol_fut_df <- gather(gfdl_seasonol_fut_df, model, prob, 3:17)
gfdl_seasonol_fut_df <- separate(data = gfdl_seasonol_fut_df, col = model, into = c("species", "season"), sep = "_gfdl_fut_")
gfdl_seasonol_fut_df$season <- factor(gfdl_seasonol_fut_df$season, levels = c("davidson", "upwelling","oceanic"))

# ggplot(gfdl_seasonol_fut_df) + geom_tile(aes(x=x,y=y, fill = prob)) + facet_wrap(~season+species, ncol = 5) +
#   geom_sf(data = coastline_sf, fill = 'black') +
#   geom_sf(data = states, fill = 'white', color = 'black') +
#   geom_sf(data = sanc, col = 'white', fill = 'transparent') +
#   geom_sf(data = ca, fill = 'transparent', color = 'white') +
#   geom_sf(data = or, fill = 'transparent', color = 'white') +
#   scale_x_continuous(limits = c(-128, -117), expand = c(0,0), breaks = c(-127, -124, -121, -118)) +
#   scale_y_continuous(limits = c(32, 47.25), expand = c(0,0)) +
#   scale_fill_viridis_c() +
#   guides(fill = 'none') +
#   # xlim(-126, -117.5) + ylim(30, 48) +  
#   labs(x = NULL, y = NULL) + theme_Publication()

# Plot IPSL ####
ipsl_seasonol_fut_stack <-raster::stack(bfal_ipsl_fut_davidson, bfal_ipsl_fut_upwelling, bfal_ipsl_fut_oceanic,
                                         caau_ipsl_fut_davidson, caau_ipsl_fut_upwelling, caau_ipsl_fut_oceanic,
                                         comu_ipsl_fut_davidson, comu_ipsl_fut_upwelling, comu_ipsl_fut_oceanic,
                                         rhau_ipsl_fut_davidson, rhau_ipsl_fut_upwelling, rhau_ipsl_fut_oceanic,
                                         sosh_ipsl_fut_davidson, sosh_ipsl_fut_upwelling, sosh_ipsl_fut_oceanic)

names(ipsl_seasonol_fut_stack) <- c('bfal_ipsl_fut_davidson', 'bfal_ipsl_fut_upwelling', 'bfal_ipsl_fut_oceanic',
                                     'caau_ipsl_fut_davidson', 'caau_ipsl_fut_upwelling', 'caau_ipsl_fut_oceanic',
                                     'comu_ipsl_fut_davidson', 'comu_ipsl_fut_upwelling', 'comu_ipsl_fut_oceanic',
                                     'rhau_ipsl_fut_davidson', 'rhau_ipsl_fut_upwelling', 'rhau_ipsl_fut_oceanic',
                                     'sosh_ipsl_fut_davidson', 'sosh_ipsl_fut_upwelling', 'sosh_ipsl_fut_oceanic')

ipsl_seasonol_fut_df <- as.data.frame(rasterToPoints(ipsl_seasonol_fut_stack, na.rm=FALSE))


ipsl_seasonol_fut_df <- gather(ipsl_seasonol_fut_df, model, prob, 3:17)
ipsl_seasonol_fut_df <- separate(data = ipsl_seasonol_fut_df, col = model, into = c("species", "season"), sep = "_ipsl_fut_")
ipsl_seasonol_fut_df$season <- factor(ipsl_seasonol_fut_df$season, levels = c("davidson", "upwelling","oceanic"))

# ggplot(ipsl_seasonol_fut_df) + geom_tile(aes(x=x,y=y, fill = prob)) + facet_wrap(~season+species, ncol = 5) +
#   geom_sf(data = coastline_sf, fill = 'black') +
#   geom_sf(data = states, fill = 'white', color = 'black') +
#   geom_sf(data = sanc, col = 'black', fill = 'transparent') +
#   geom_sf(data = ca, fill = 'transparent', color = 'black') +
#   geom_sf(data = or, fill = 'transparent', color = 'black') +
#   scale_x_continuous(limits = c(-128, -117), expand = c(0,0), breaks = c(-127, -124, -121, -118)) +
#   scale_y_continuous(limits = c(32, 47.25), expand = c(0,0)) +
#   scale_fill_viridis_c() +
#   guides(fill = 'none') +
#   # xlim(-126, -117.5) + ylim(30, 48) +  
#   labs(x = NULL, y = NULL) + theme_Publication()

# Plot HAD ####
had_seasonol_fut_stack <-raster::stack(bfal_had_fut_davidson, bfal_had_fut_upwelling, bfal_had_fut_oceanic,
                                        caau_had_fut_davidson, caau_had_fut_upwelling, caau_had_fut_oceanic,
                                        comu_had_fut_davidson, comu_had_fut_upwelling, comu_had_fut_oceanic,
                                        rhau_had_fut_davidson, rhau_had_fut_upwelling, rhau_had_fut_oceanic,
                                        sosh_had_fut_davidson, sosh_had_fut_upwelling, sosh_had_fut_oceanic)

names(had_seasonol_fut_stack) <- c('bfal_had_fut_davidson', 'bfal_had_fut_upwelling', 'bfal_had_fut_oceanic',
                                    'caau_had_fut_davidson', 'caau_had_fut_upwelling', 'caau_had_fut_oceanic',
                                    'comu_had_fut_davidson', 'comu_had_fut_upwelling', 'comu_had_fut_oceanic',
                                    'rhau_had_fut_davidson', 'rhau_had_fut_upwelling', 'rhau_had_fut_oceanic',
                                    'sosh_had_fut_davidson', 'sosh_had_fut_upwelling', 'sosh_had_fut_oceanic')

had_seasonol_fut_df <- as.data.frame(rasterToPoints(had_seasonol_fut_stack, na.rm=FALSE))


had_seasonol_fut_df <- gather(had_seasonol_fut_df, model, prob, 3:17)
had_seasonol_fut_df <- separate(data = had_seasonol_fut_df, col = model, into = c("species", "season"), sep = "_had_fut_")
had_seasonol_fut_df$season <- factor(had_seasonol_fut_df$season, levels = c("davidson", "upwelling","oceanic"))

# ggplot(had_seasonol_fut_df) + geom_tile(aes(x=x,y=y, fill = prob)) + facet_wrap(~season+species, ncol = 5) +
#   geom_sf(data = coastline_sf, fill = 'black') +
#   geom_sf(data = states, fill = 'white', color = 'black') +
#   geom_sf(data = sanc, col = 'black', fill = 'transparent') +
#   geom_sf(data = ca, fill = 'transparent', color = 'black') +
#   geom_sf(data = or, fill = 'transparent', color = 'black') +
#   scale_x_continuous(limits = c(-128, -117), expand = c(0,0), breaks = c(-127, -124, -121, -118)) +
#   scale_y_continuous(limits = c(32, 47.25), expand = c(0,0)) +
#   scale_fill_viridis_c() +
#   guides(fill = 'none') +
#   # xlim(-126, -117.5) + ylim(30, 48) +  
#   labs(x = NULL, y = NULL) + theme_Publication()


# Plot Ensemble ####
gfdl_seasonol_fut_df$mod <- 'gfdl'
ipsl_seasonol_fut_df$mod <- 'ipsl'
had_seasonol_fut_df$mod <- 'had'

ensem_fut_seasonclimat_df <- cbind(gfdl_seasonol_fut_df, ipsl_seasonol_fut_df[,c(5,6)], had_seasonol_fut_df[,c(5,6)])
ensem_fut_seasonclimat_df <- ensem_fut_seasonclimat_df[!duplicated(as.list(ensem_fut_seasonclimat_df))]

ensem_fut_seasonclimat_df <- ensem_fut_seasonclimat_df[,-c(6, 8, 10)]
head(ensem_fut_seasonclimat_df)
names(ensem_fut_seasonclimat_df)[5:7] <- c('prob_gfdl', 'prob_ipsl', 'prob_had')

ensem_fut_seasonclimat_df$prob_em <- rowMeans(ensem_fut_seasonclimat_df[5:7], na.rm=TRUE) # calc ensemble mean
ensem_fut_seasonclimat_df$prob_esd <- apply(ensem_fut_seasonclimat_df[,5:7],1,sd)  # calc ensemble sd
head(ensem_fut_seasonclimat_df)
ensem_fut_seasonclimat_df$dist <- raster::extract(dist, ensem_fut_seasonclimat_df[,1:2])

# ensemble mean future HSI
ggplot(subset(ensem_fut_seasonclimat_df, dist < 273000)) + 
  geom_tile(aes(x=x,y=y, fill = prob_em)) + 
  facet_wrap(~season+species, ncol = 5) +
  geom_sf(data = coastline_sf, fill = 'black') +
  geom_sf(data = states, fill = 'white', color = 'black') +
  geom_sf(data = sanc, col = 'white', fill = 'transparent') +
  geom_sf(data = ca, fill = 'transparent', color = 'white') +
  geom_sf(data = or, fill = 'transparent', color = 'white') +
  scale_x_continuous(limits = c(-128, -117), expand = c(0,0), breaks = c(-127, -124, -121, -118)) +
  scale_y_continuous(limits = c(32, 47.25), expand = c(0,0)) +
  scale_fill_viridis_c() +
  # guides(fill = 'none') +
  labs(x = NULL, y = NULL, fill = 'HSI') + 
  theme_Publication() + theme(legend.position = c(0.95,0.95), legend.direction = 'vertical')

# ensemble sd future HSI
ggplot(subset(ensem_fut_seasonclimat_df, dist < 273000)) + 
  geom_tile(aes(x=x,y=y, fill = prob_esd)) + 
  facet_wrap(~season+species, ncol = 5) +
  geom_sf(data = coastline_sf, fill = 'black') +
  geom_sf(data = states, fill = 'white', color = 'black') +
  geom_sf(data = sanc, col = 'white', fill = 'transparent') +
  geom_sf(data = ca, fill = 'transparent', color = 'white') +
  geom_sf(data = or, fill = 'transparent', color = 'white') +
  scale_x_continuous(limits = c(-128, -117), expand = c(0,0), breaks = c(-127, -124, -121, -118)) +
  scale_y_continuous(limits = c(32, 47.25), expand = c(0,0)) +
  scale_fill_viridis_c() +
  # guides(fill = 'none') +
  labs(x = NULL, y = NULL, fill = 'ESM SD') + 
  theme_Publication() + theme(legend.position = c(0.95,0.95), legend.direction = 'vertical')

#_____ ####
# COMBINE & PLOT ####
ensem_fut_seasonclimat_df$period <- 'future'
ensem_hist_seasonclimat_df$period <- 'historical'

ensem_hist_fut_seasonclimat_df <- cbind(ensem_hist_seasonclimat_df, ensem_fut_seasonclimat_df)
ensem_hist_fut_seasonclimat_df <- ensem_hist_fut_seasonclimat_df[!duplicated(as.list(ensem_hist_fut_seasonclimat_df))]
head(ensem_hist_fut_seasonclimat_df)
names(ensem_hist_fut_seasonclimat_df)  <- c("x", "y", "species", "season", 
                                            "prob_gfdl_hist",   "prob_ipsl_hist" ,  "prob_had_hist" ,   "prob_em_hist",     "prob_esd_hist" ,   
                                            "dist" , 'period', 
                                            "prob_gfdl_fut", "prob_ipsl_fut", "prob_had_fut", "prob_em_fut" , "prob_esd_fut", 
                                            "period.1" )

ensem_hist_fut_seasonclimat_df <- ensem_hist_fut_seasonclimat_df[,-c(11,17)]
head(ensem_hist_fut_seasonclimat_df)

ensem_hist_fut_seasonclimat_df$perc_change_gfdl <- (((ensem_hist_fut_seasonclimat_df$prob_gfdl_fut - ensem_hist_fut_seasonclimat_df$prob_gfdl_hist)/ensem_hist_fut_seasonclimat_df$prob_gfdl_hist)*100)
ensem_hist_fut_seasonclimat_df$perc_change_ipsl <- (((ensem_hist_fut_seasonclimat_df$prob_ipsl_fut - ensem_hist_fut_seasonclimat_df$prob_ipsl_hist)/ensem_hist_fut_seasonclimat_df$prob_ipsl_hist)*100)
ensem_hist_fut_seasonclimat_df$perc_change_had <- (((ensem_hist_fut_seasonclimat_df$prob_had_fut - ensem_hist_fut_seasonclimat_df$prob_had_hist)/ensem_hist_fut_seasonclimat_df$prob_had_hist)*100)
names(ensem_hist_fut_seasonclimat_df)
ensem_hist_fut_seasonclimat_df$perc_change_em <- rowMeans(ensem_hist_fut_seasonclimat_df[14:16], na.rm=TRUE) # calc ensemble mean
ensem_hist_fut_seasonclimat_df$perc_change_esd <- apply(ensem_hist_fut_seasonclimat_df[,14:16],1,sd)  # calc ensemble sd

ensem_hist_fut_seasonclimat_df <- ensem_hist_fut_seasonclimat_df %>%
  group_by(species) %>% # normalize by species
  mutate(norm_prob_em_hist = (prob_em_hist - min(prob_em_hist)) / (max(prob_em_hist) - min(prob_em_hist))) %>%
  ungroup()  

ensem_hist_fut_seasonclimat_df$species <- factor(ensem_hist_fut_seasonclimat_df$species, 
                              levels = c('bfal', 'sosh', 'comu', 'caau', 'rhau'),
                              labels = c('BFAL', 'SOSH', 'COMU', 'CAAU', 'RHAU'))

# Check the changes
levels(ensem_hist_fut_seasonclimat_df$species)
unique(ensem_hist_fut_seasonclimat_df$species)

# Plot percent change in ensemble mean HSI from historical to future in each season * species
p4 <- ggplot(subset(ensem_hist_fut_seasonclimat_df, dist < 273000)) + 
  geom_tile(aes(x=x,y=y, fill = perc_change_em, alpha = norm_prob_em_hist)) + 
  facet_wrap(~season+species, ncol = 5) +
  geom_sf(data = coastline_sf) +
  geom_sf(data = states, fill = 'white', color = 'black') +
  geom_sf(data = sanc, col = 'black', fill = 'transparent') +
  geom_sf(data = chnms, col = 'black', fill = 'transparent') +
  geom_sf(data = ca, color = 'black', fill = 'white', alpha = 0.9) +
  geom_sf(data = or, color = 'black', fill = 'white', alpha = 0.9) +
  scale_x_continuous(limits = c(-128, -117), expand = c(0,0), breaks = c(-126, -122, -118)) +
  scale_y_continuous(limits = c(32, 47.25), expand = c(0,0)) +
  scale_fill_gradient2(low = 'red', mid = 'grey', high = 'blue', midpoint = 0,
                       breaks = c(-15, -7.5, 0, 7.5)) + 
  scale_alpha(range = c(0, 1), guide = 'legend', name = 'Historical\nHSI') + 
  labs(x = NULL, y = NULL, fill = '%\nChange\nHSI', alpha = 'Historical\nHSI') + 
  theme_Publication() + theme(legend.position = 'right',
                              legend.box = 'vertical', 
                              panel.border = element_rect(color = "black", fill = NA),  
                              panel.spacing = unit(0, "cm"),
                              legend.direction = 'vertical')

p4 

# Plot percent change in ensemble sd HSI from historical to future in each season * species
ggplot(subset(ensem_hist_fut_seasonclimat_df, dist < 273000)) + 
  geom_tile(aes(x=x,y=y, fill = perc_change_esd)) +
  facet_wrap(~season+species, ncol = 5) +
  geom_sf(data = coastline_sf, fill = 'black') +
  geom_sf(data = states, fill = 'white', color = 'black') +
  geom_sf(data = sanc, col = 'white', alpha = 0) +
  geom_sf(data = ca, color = 'white',  alpha = 0) +
  geom_sf(data = or, color = 'white',  alpha = 0) +
  scale_x_continuous(limits = c(-128, -117), expand = c(0,0), breaks = c(-126, -122, -118)) +
  scale_y_continuous(limits = c(32, 47.25), expand = c(0,0)) +
  scale_fill_viridis_c() +  labs(x = NULL, y = NULL, fill = '%\nChange\nSD') + 
  theme_Publication() + theme(legend.position  = c(0.95,0.93),
                              legend.box = 'vertical', 
                              panel.border = element_rect(color = "black", fill = NA),  
                              panel.spacing = unit(0, "cm"),
                              legend.direction = 'vertical')
