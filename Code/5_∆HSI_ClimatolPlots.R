# HEADER ####
# Author: Ryan Gasbarro; rygasbar@ucsc.edu
# Script Purpose: Plot historical ensemble mean and change in HSI between future (2070-2100) & historical (1985-2015) time periods for each seabird spp

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
library(ggrepel)
library(gdata)
library(viridis)
library(viridisLite)
library(biscale)
library(sf)
library(geosphere)
library(rmapshaper)
library(rnaturalearth)
library(rnaturalearthdata)
library(wesanderson)

#optional theme for plotting
PubTheme <- function(base_size=12, base_family="helvetica") {
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
            axis.text = element_text(size = 10), 
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

# Load coastline, sanctuaries, weas, etc. shapefiles ####
data.dir <- "~/SeabirdModels/Data" # * CHANGE TO YOUR DATA DIRECTORY
setwd(data.dir)
z <- as.data.frame(rasterToPoints(raster('z.grd')))
z_rast <- raster('z.grd')
z_rast <- round(z_rast)
z$z2 <- round(z$layer)
names(z) <- c('lon', 'lat', 'z','z_rounded')
head(z)

dist <- raster('Distance_from_shore.grd')
dist_df <- as.data.frame(rasterToPoints(raster('Distance_from_shore.grd')))

coastline <- ne_countries(scale = "medium", returnclass = 'sf')
coastline_bbox <- c(xmin=-134,xmax=-116,ymin=30,ymax=50)
coastline_sf <- coastline %>% 
  filter(admin %in% c("Canada", "Mexico", "United States of America")) %>% 
  st_crop(coastline_bbox)
states <- ne_states(country = "united states of america", returnclass = "sf")
states <- states %>%
  filter(name %in% c('Washington', 'Oregon', 'California', 'Nevada', 'Idaho')) %>% 
  st_crop(coastline_bbox)


shape.dir <- '~/SeabirdModels/Shapefiles' # directory where National Marine Sanctuary (NMS) shapefiles are stored
setwd(shape.dir)

sanc <- read_sf('Combined')
sanc$attribute1[4] <- 'CBNMS' 
chnms <- read_sf('CHNMS_Boundary_10152024')

wind.dir <- '~/SeabirdModels/Shapefiles/WindAreas' # directory where wind area shapefiles are stored
setwd(wind.dir)

ca <- read_sf('CA_FinalLeaseAreas_Shapefiles/LeaseAreas_CA_Outlines_2022_09_08.shp')
or <- read_sf('OR_WEA_GIS_Data/CDR_CDR_OREGON_WEA_2024_BOEM.shp')

# Load Climatologies ####

# Load Rasters covering historical period (1985-2015) #
# load gfdl climatologies
bfal_gfdl_hist <- raster("~/SeabirdModels/ProjectionClimatologies/bfal/gfdl/Climatol_1985-2015/BFAL_OccProbClimatol_GFDL.tif")
caau_gfdl_hist <- raster("~/SeabirdModels/ProjectionClimatologies/caau/gfdl/Climatol_1985-2015/CAAU_OccProbClimatol_GFDL.tif")
comu_gfdl_hist <- raster("~/SeabirdModels/ProjectionClimatologies/comu/gfdl/Climatol_1985-2015/COMU_OccProbClimatol_GFDL.tif")
rhau_gfdl_hist <- raster("~/SeabirdModels/ProjectionClimatologies/rhau/gfdl/Climatol_1985-2015/RHAU_OccProbClimatol_GFDL.tif")
sosh_gfdl_hist <- raster("~/SeabirdModels/ProjectionClimatologies/sosh/gfdl/Climatol_1985-2015/SOSH_OccProbClimatol_GFDL.tif")

# load ipsl climatologies
bfal_ipsl_hist <- raster("~/SeabirdModels/ProjectionClimatologies/bfal/ipsl/Climatol_1985-2015/BFAL_OccProbClimatol_IPSL.tif")
caau_ipsl_hist <- raster("~/SeabirdModels/ProjectionClimatologies/caau/ipsl/Climatol_1985-2015/CAAU_OccProbClimatol_IPSL.tif")
comu_ipsl_hist <- raster("~/SeabirdModels/ProjectionClimatologies/comu/ipsl/Climatol_1985-2015/COMU_OccProbClimatol_IPSL.tif")
rhau_ipsl_hist <- raster("~/SeabirdModels/ProjectionClimatologies/rhau/ipsl/Climatol_1985-2015/RHAU_OccProbClimatol_IPSL.tif")
sosh_ipsl_hist <- raster("~/SeabirdModels/ProjectionClimatologies/sosh/ipsl/Climatol_1985-2015/SOSH_OccProbClimatol_IPSL.tif")

# load had climatologies
bfal_had_hist <- raster("~/SeabirdModels/ProjectionClimatologies/bfal/had/Climatol_1985-2015/BFAL_OccProbClimatol_HAD.tif")
caau_had_hist <- raster("~/SeabirdModels/ProjectionClimatologies/caau/had/Climatol_1985-2015/CAAU_OccProbClimatol_HAD.tif")
comu_had_hist <- raster("~/SeabirdModels/ProjectionClimatologies/comu/had/Climatol_1985-2015/COMU_OccProbClimatol_HAD.tif")
rhau_had_hist <- raster("~/SeabirdModels/ProjectionClimatologies/rhau/had/Climatol_1985-2015/RHAU_OccProbClimatol_HAD.tif")
sosh_had_hist <- raster("~/SeabirdModels/ProjectionClimatologies/sosh/had/Climatol_1985-2015/SOSH_OccProbClimatol_HAD.tif")

bfal_ensem_mean <- mean(raster::stack(bfal_gfdl_hist, bfal_ipsl_hist, bfal_had_hist)) # ens mean
caau_ensem_mean <- mean(raster::stack(caau_gfdl_hist, caau_ipsl_hist, caau_had_hist))
comu_ensem_mean <- mean(raster::stack(comu_gfdl_hist, comu_ipsl_hist, comu_had_hist))
rhau_ensem_mean <- mean(raster::stack(rhau_gfdl_hist, rhau_ipsl_hist, rhau_had_hist))
sosh_ensem_mean <- mean(raster::stack(sosh_gfdl_hist, sosh_ipsl_hist, sosh_had_hist))

# Load Rasters Future (2070-2100) #

# load gfdl climatologies
bfal_gfdl_fut <- raster("~/SeabirdModels/ProjectionClimatologies/bfal/gfdl/Climatol_2070-2100/BFAL_OccProbClimatol_GFDL.tif")
caau_gfdl_fut <- raster("~/SeabirdModels/ProjectionClimatologies/caau/gfdl/Climatol_2070-2100/CAAU_OccProbClimatol_GFDL.tif")
comu_gfdl_fut <- raster("~/SeabirdModels/ProjectionClimatologies/comu/gfdl/Climatol_2070-2100/COMU_OccProbClimatol_GFDL.tif")
rhau_gfdl_fut <- raster("~/SeabirdModels/ProjectionClimatologies/rhau/gfdl/Climatol_2070-2100/RHAU_OccProbClimatol_GFDL.tif")
sosh_gfdl_fut <- raster("~/SeabirdModels/ProjectionClimatologies/sosh/gfdl/Climatol_2070-2100/SOSH_OccProbClimatol_GFDL.tif")

# load ipsl climatologies
bfal_ipsl_fut <- raster("~/SeabirdModels/ProjectionClimatologies/bfal/ipsl/Climatol_2070-2100/BFAL_OccProbClimatol_IPSL.tif")
caau_ipsl_fut <- raster("~/SeabirdModels/ProjectionClimatologies/caau/ipsl/Climatol_2070-2100/CAAU_OccProbClimatol_IPSL.tif")
comu_ipsl_fut <- raster("~/SeabirdModels/ProjectionClimatologies/comu/ipsl/Climatol_2070-2100/COMU_OccProbClimatol_IPSL.tif")
rhau_ipsl_fut <- raster("~/SeabirdModels/ProjectionClimatologies/rhau/ipsl/Climatol_2070-2100/RHAU_OccProbClimatol_IPSL.tif")
sosh_ipsl_fut <- raster("~/SeabirdModels/ProjectionClimatologies/sosh/ipsl/Climatol_2070-2100/SOSH_OccProbClimatol_IPSL.tif")

# load had climatologies
bfal_had_fut <- raster("~/SeabirdModels/ProjectionClimatologies/bfal/had/Climatol_2070-2100/BFAL_OccProbClimatol_HAD.tif")
caau_had_fut <- raster("~/SeabirdModels/ProjectionClimatologies/caau/had/Climatol_2070-2100/CAAU_OccProbClimatol_HAD.tif")
comu_had_fut <- raster("~/SeabirdModels/ProjectionClimatologies/comu/had/Climatol_2070-2100/COMU_OccProbClimatol_HAD.tif")
rhau_had_fut <- raster("~/SeabirdModels/ProjectionClimatologies/rhau/had/Climatol_2070-2100/RHAU_OccProbClimatol_HAD.tif")
sosh_had_fut <- raster("~/SeabirdModels/ProjectionClimatologies/sosh/had/Climatol_2070-2100/SOSH_OccProbClimatol_HAD.tif")

bfal_ensem_mean_fut <- mean(raster::stack(bfal_gfdl_fut, bfal_ipsl_fut, bfal_had_fut))
caau_ensem_mean_fut <- mean(raster::stack(caau_gfdl_fut, caau_ipsl_fut, caau_had_fut))
comu_ensem_mean_fut <- mean(raster::stack(comu_gfdl_fut, comu_ipsl_fut, comu_had_fut))
rhau_ensem_mean_fut <- mean(raster::stack(rhau_gfdl_fut, rhau_ipsl_fut, rhau_had_fut))
sosh_ensem_mean_fut <- mean(raster::stack(sosh_gfdl_fut, sosh_ipsl_fut, sosh_had_fut))

# calculate changes in ensemble mean HSI between time periods:
bfal_ensem_mean_chng <- bfal_ensem_mean_fut - bfal_ensem_mean
caau_ensem_mean_chng <- caau_ensem_mean_fut - caau_ensem_mean
comu_ensem_mean_chng <- comu_ensem_mean_fut - comu_ensem_mean
rhau_ensem_mean_chng <- rhau_ensem_mean_fut - rhau_ensem_mean
sosh_ensem_mean_chng <- sosh_ensem_mean_fut - sosh_ensem_mean


hist_stack <-             raster::stack(bfal_ipsl_hist, bfal_gfdl_hist, bfal_had_hist,
                                        caau_ipsl_hist, caau_gfdl_hist, caau_had_hist,
                                        comu_ipsl_hist, comu_gfdl_hist, comu_had_hist,
                                        rhau_ipsl_hist, rhau_gfdl_hist, rhau_had_hist,
                                        sosh_ipsl_hist, sosh_gfdl_hist, sosh_had_hist)

names(hist_stack) <-             c('bfal_ipsl_hist', 'bfal_gfdl_hist', 'bfal_had_hist',
                                    'caau_ipsl_hist', 'caau_gfdl_hist', 'caau_had_hist',
                                    'comu_ipsl_hist', 'comu_gfdl_hist', 'comu_had_hist',
                                    'rhau_ipsl_hist', 'rhau_gfdl_hist', 'rhau_had_hist',
                                    'sosh_ipsl_hist', 'sosh_gfdl_hist', 'sosh_had_hist')

hist_df <- as.data.frame(rasterToPoints(hist_stack, na.rm=FALSE))
hist_df <- gather(hist_df, model, prob, 3:17)
hist_df <- separate(data = hist_df, col = model, into = c("species", "mod",'period'), sep = "_")
head(hist_df)

fut_stack <-             raster::stack(bfal_ipsl_fut, bfal_gfdl_fut, bfal_had_fut,
                                        caau_ipsl_fut, caau_gfdl_fut, caau_had_fut,
                                        comu_ipsl_fut, comu_gfdl_fut, comu_had_fut,
                                        rhau_ipsl_fut, rhau_gfdl_fut, rhau_had_fut,
                                        sosh_ipsl_fut, sosh_gfdl_fut, sosh_had_fut)

names(fut_stack) <-             c('bfal_ipsl_fut', 'bfal_gfdl_fut', 'bfal_had_fut',
                                   'caau_ipsl_fut', 'caau_gfdl_fut', 'caau_had_fut',
                                   'comu_ipsl_fut', 'comu_gfdl_fut', 'comu_had_fut',
                                   'rhau_ipsl_fut', 'rhau_gfdl_fut', 'rhau_had_fut',
                                   'sosh_ipsl_fut', 'sosh_gfdl_fut', 'sosh_had_fut')

fut_df <- as.data.frame(rasterToPoints(fut_stack, na.rm=FALSE))
fut_df <- gather(fut_df, model, prob, 3:17)
fut_df <- separate(data = fut_df, col = model, into = c("species", "mod",'period'), sep = "_")
head(fut_df)

delta_stack <-             raster::stack(bfal_ensem_mean_chng, 
                                         caau_ensem_mean_chng, 
                                         comu_ensem_mean_chng, 
                                         rhau_ensem_mean_chng, 
                                         sosh_ensem_mean_chng) 

names(delta_stack) <-             c('bfal_ensem_delta', 
                                  'caau_ensem_delta', 
                                  'comu_ensem_delta', 
                                  'rhau_ensem_delta',
                                  'sosh_ensem_delta')

delta_df <- as.data.frame(rasterToPoints(delta_stack, na.rm=FALSE))
delta_df <- gather(delta_df, model, prob, 3:7)
delta_df <- separate(data = delta_df, col = model, into = c("species", "mod",'period'), sep = "_")
head(delta_df)

ensem_hist_stack <-             raster::stack(bfal_ensem_mean, 
                                              caau_ensem_mean, 
                                              comu_ensem_mean, 
                                              rhau_ensem_mean, 
                                              sosh_ensem_mean) 

names(ensem_hist_stack) <-          c('bfal_ensem_hist', 
                                    'caau_ensem_hist', 
                                    'comu_ensem_hist', 
                                    'rhau_ensem_hist',
                                    'sosh_ensem_hist')

ensem_hist_df <- as.data.frame(rasterToPoints(ensem_hist_stack, na.rm=FALSE))
ensem_hist_df <- gather(ensem_hist_df, model, prob, 3:7)
ensem_hist_df <- separate(data = ensem_hist_df, col = model, into = c("species", "mod",'period'), sep = "_")
head(ensem_hist_df)

df <- rbind(hist_df, fut_df, delta_df, ensem_hist_df) # combine all dataframes

df$dist <- raster::extract(dist, df[,1:2]) # distance from shore
head(df) # manuall inspect

# PLOT ####
pal <- wes_palette("Zissou1", 100, type = 'continuous')

df$species <- factor(df$species, 
                             levels = c('bfal', 'sosh', 'comu', 'caau', 'rhau'),
                             labels = c('BFAL', 'SOSH', 'COMU', 'CAAU', 'RHAU'))


a <- ggplot() + 
  geom_tile(data = subset(df, period == 'hist' & mod == 'ensem' & dist < 273000), aes(x = x, y = y, fill = prob)) +
  facet_wrap(~species, ncol = 5) +
  geom_sf(data = coastline_sf, fill = 'transparent', col = 'black') + 
  geom_sf(data = states, fill = 'transparent', color = 'black') +
  geom_sf(data = sanc, col = 'black', fill = 'transparent', lwd = 0.52) +
  geom_sf(data = chnms, col = 'black', fill = 'transparent', lwd = 0.52) +
  geom_sf(data = ca, fill = 'white', color = 'black', lwd = 0.52, alpha = 0.8) +
  geom_sf(data = or, fill = 'white', color = 'black', lwd = 0.52, alpha = 0.95) +
  scale_x_continuous(limits = c(-128, -117), expand = c(0,0), breaks = c(-128, -125, -122,  -119)) +
  scale_y_continuous(limits = c(32, 47.25), expand = c(0,0)) +
  labs(x = NULL, y = NULL, fill = 'Habitat \nSuitability \nIndex [HSI]', shape = 'Platform') +
  PubTheme() + 
  scale_fill_gradientn(colours = pal) +
  theme(legend.position = 'inside',
        legend.position.inside = c(0.15, 0.75), 
        legend.background = element_rect(fill = 'white', colour = 1),
        legend.direction = 'vertical',
        strip.background = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = margin(2, 1, 0, 2),
        strip.text.x = element_blank())

b <- ggplot() + 
  geom_tile(data = subset(df, period == 'delta' & dist < 273000), aes(x = x, y = y, fill = prob)) +
  facet_wrap(~species, ncol = 5) +
  geom_sf(data = coastline_sf, fill = 'transparent', col = 'black') + 
  geom_sf(data = states, fill = 'transparent', color = 'black') +
  geom_sf(data = sanc, col = 'black', fill = 'transparent', lwd = 0.52) +
  geom_sf(data = chnms, col = 'black', fill = 'transparent', lwd = 0.52) +
  geom_sf(data = ca, fill = 'white', color = 'black', lwd = 0.52, alpha = 0.8) +
  geom_sf(data = or, fill = 'white', color = 'black', lwd = 0.52, alpha = 0.95) +
  scale_x_continuous(limits = c(-128, -117), expand = c(0,0), breaks = c(-128,-125, -122,  -119)) +
  scale_y_continuous(limits = c(32, 47.25), expand = c(0,0)) +
  labs(x = NULL, y = NULL, fill = '∆HSI', shape = 'Platform') +
  PubTheme() + 
  scale_fill_gradient2(low = "red", mid = "grey",  high = "blue",  midpoint = 0) +
  theme(legend.position = 'inside',
        legend.position.inside = c(0.15, 0.75), 
        legend.background = element_rect(fill = 'white', colour = 1),
        legend.direction = 'vertical',
        strip.background = element_blank(),
        axis.text.y = element_blank(),
        plot.margin = margin(2, 1, 0, 2),
        strip.text.x = element_blank())

grid.arrange(a,b, padding = unit(0.25, "mm"))

