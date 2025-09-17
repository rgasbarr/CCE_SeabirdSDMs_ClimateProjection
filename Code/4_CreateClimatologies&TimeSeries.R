# HEADER ####
# Author: Ryan Gasbarro; rygasbar@ucsc.edu
# Purpose: Example script detailing process of creating habitat suitability climatologies & time-series 
# This example uses one species: black-footed albatross (BFAL)
# And one earth-system model (IPSL)
# See Gasbarro et al. (2025) PLoS Climate for more information

# _____________________#####
# CREATE CLIMATOLOGIES ####
# 1985-2015 IPSL ####
predDir="bfal/dailyP_IPSL/"  # change to the folder where your daily projections are stored
setwd(predDir)
fpath <- getwd()

# List all files in the directory
raster_files <- list.files(pattern = "\\.gri$")

extract_date <- function(filename) {
  date_string <- gsub("bfal_", "", filename) # Remove 'bfal_' prefix
  date_string <- gsub("_mean\\.gri", "", date_string) # Remove '_mean.gri' suffix
  return(ymd(date_string)) # Parse date
}

# Filter 1985 to 2015
raster_files <- raster_files[extract_date(raster_files) >= ymd("1985-01-01") & 
                               extract_date(raster_files) <= ymd("2015-12-31")]

means_folder <- file.path(fpath, "Climatol_1985-2015")
if (!dir.exists(means_folder)) {
  dir.create(means_folder)
}

rasters <- lapply(raster_files, raster)
stack_data <- raster::stack(rasters)
mean_r <- mean(stack_data)

writeRaster(mean_r, 
            filename = "Climatol_1985-2015/BFAL_OccProbClimatol_IPSL.grd", 
            format = "raster", overwrite = TRUE)

rm(list = ls()) 
# 2035-2065 IPSL ####
predDir="bfal/dailyP_IPSL/"
setwd(predDir)
fpath <- getwd()

# List all files in the directory
raster_files <- list.files(pattern = "\\.gri$")

extract_date <- function(filename) {
  date_string <- gsub("bfal_", "", filename) # Remove 'bfal_' prefix
  date_string <- gsub("_mean\\.gri", "", date_string) # Remove '_mean.gri' suffix
  return(ymd(date_string)) # Parse date
}

# Filter 2035 to 2065
raster_files <- raster_files[extract_date(raster_files) >= ymd("2035-01-01") & 
                               extract_date(raster_files) <= ymd("2065-12-31")]

means_folder <- file.path(fpath, "Climatol_2035-2065")
if (!dir.exists(means_folder)) {
  dir.create(means_folder)
}

rasters <- lapply(raster_files, raster)
stack_data <- raster::stack(rasters)
mean_r <- mean(stack_data)
my.at <- seq(0, 1, 0.01)

writeRaster(mean_r, 
            filename = "Climatol_2035-2065/BFAL_OccProbClimatol_IPSL.grd", 
            format = "raster", overwrite = TRUE)

rm(list = ls())

# 2070-2100 IPSL ####
predDir="bfal/dailyP_IPSL/"  # change to the folder where your daily projections are stored
setwd(predDir)
fpath <- getwd()

# List all files in the directory
raster_files <- list.files(pattern = "\\.gri$")

extract_date <- function(filename) {
  date_string <- gsub("bfal_", "", filename) # Remove 'bfal_' prefix
  date_string <- gsub("_mean\\.gri", "", date_string) # Remove '_mean.gri' suffix
  return(ymd(date_string)) # Parse date
}

# Filter 2070 to 2100
raster_files <- raster_files[extract_date(raster_files) >= ymd("2070-01-01") & 
                               extract_date(raster_files) <= ymd("2100-12-31")]

means_folder <- file.path(fpath, "Climatol_2070-2100")
if (!dir.exists(means_folder)) {
  dir.create(means_folder)
}

rasters <- lapply(raster_files, raster)
stack_data <- raster::stack(rasters)
mean_r <- mean(stack_data)
my.at <- seq(0, 1, 0.01)

writeRaster(mean_r, 
            filename = "Climatol_2070-2100/BFAL_OccProbClimatol_IPSL.grd", 
            format = "raster", overwrite = TRUE)

rm(list = ls())

# _____________________#####
# ENSEMBLE CLIMATOLOGIES ####
# Ensemble Mean & SD 1985-2015 ####
# NOTE YOU WILL HAVE TO DO THE ABOVE PROCESS OF CREATING CLIMATOLOGIES
# FOR THE OTHER TWO ESMs (i.e. GFDL, HAD) before proceeding to calculate ensemble means) !!!!!!!!!!!!!!!

ipsl_dir <- predDir
gfdl_dir <- "bfal/dailyP_GFDL/"
had_dir <- "bfal/dailyP_HAD/"

setwd(gfdl_dir)
bfal_gfdl_mean_hist <- raster("Climatol_1985-2015/BFAL_OccProbClimatol_GFDL.tif")

setwd(ipsl_dir)
bfal_ipsl_mean_hist <- raster("Climatol_1985-2015/BFAL_OccProbClimatol_IPSL.tif")

setwd(had_dir)
bfal_had_mean_hist <- raster("Climatol_1985-2015/BFAL_OccProbClimatol_HAD.tif")

bfal_ESMstack_hist <- raster::stack(bfal_gfdl_mean_hist, bfal_ipsl_mean_hist, bfal_had_mean_hist)

bfal_ESMmean_hist <- mean(bfal_ESMstack_hist)
bfal_ESMsd_hist <- calc(bfal_ESMstack_hist, fun = sd)

ensem.dir <- 'bfal/Ensemble' # changed to where you would like ensemble mean & st. dev. files to be stored
setwd(ensem.dir)
writeRaster(bfal_ESMmean_hist, filename = "BFAL_ESMmean_1985-2015.tif", format = "GTiff", overwrite = TRUE)
writeRaster(bfal_ESMsd_hist, filename = "BFAL_ESMsd_1985-2015.tif", format = "GTiff", overwrite = TRUE)

# Ensemble Mean & SD 2035-2065 ####
# NOTE YOU WILL HAVE TO DO THE ABOVE PROCESS OF CREATING CLIMATOLOGIES
# FOR THE OTHER TWO ESMs (i.e. GFDL, HAD) before proceeding to calculate ensemble means) !!!!!!!!!!!!!!!
ipsl_dir <- predDir
gfdl_dir <- "bfal/dailyP_GFDL/"
had_dir <- "bfal/dailyP_HAD/"

setwd(gfdl_dir)
bfal_gfdl_mean_hist <- raster("Climatol_2035-2065/BFAL_OccProbClimatol_GFDL.tif")

setwd(ipsl_dir)
bfal_ipsl_mean_hist <- raster("Climatol_2035-2065/BFAL_OccProbClimatol_IPSL.tif")

setwd(had_dir)
bfal_had_mean_hist <- raster("Climatol_2035-2065/BFAL_OccProbClimatol_HAD.tif")

bfal_ESMstack_hist <- raster::stack(bfal_gfdl_mean_hist, bfal_ipsl_mean_hist, bfal_had_mean_hist)

bfal_ESMmean_hist <- mean(bfal_ESMstack_hist)
bfal_ESMsd_hist <- calc(bfal_ESMstack_hist, fun = sd)

ensem.dir <- 'bfal/Ensemble' # changed to where you would like ensemble mean & st. dev. files to be stored
setwd(ensem.dir)
writeRaster(bfal_ESMmean_hist, filename = "BFAL_ESMmean_2035-2065.tif", format = "GTiff", overwrite = TRUE)
writeRaster(bfal_ESMsd_hist, filename = "BFAL_ESMsd_2035-2065.tif", format = "GTiff", overwrite = TRUE)

# Ensemble Mean & SD 2070-2100 ####
# NOTE YOU WILL HAVE TO DO THE ABOVE PROCESS OF CREATING CLIMATOLOGIES
# FOR THE OTHER TWO ESMs (i.e. GFDL, HAD) before proceeding to calculate ensemble means) !!!!!!!!!!!!!!!

ipsl_dir <- predDir
gfdl_dir <- "bfal/dailyP_GFDL/"
had_dir <- "bfal/dailyP_HAD/"

setwd(gfdl_dir)
bfal_gfdl_mean_hist <- raster("Climatol_2070-2100/BFAL_OccProbClimatol_GFDL.tif")

setwd(ipsl_dir)
bfal_ipsl_mean_hist <- raster("Climatol_2070-2100/BFAL_OccProbClimatol_IPSL.tif")

setwd(had_dir)
bfal_had_mean_hist <- raster("Climatol_2070-2100/BFAL_OccProbClimatol_HAD.tif")

bfal_ESMstack_hist <- raster::stack(bfal_gfdl_mean_hist, bfal_ipsl_mean_hist, bfal_had_mean_hist)

bfal_ESMmean_hist <- mean(bfal_ESMstack_hist)
bfal_ESMsd_hist <- calc(bfal_ESMstack_hist, fun = sd)

ensem.dir <- 'bfal/Ensemble' # changed to where you would like ensemble mean & st. dev. files to be stored
setwd(ensem.dir)
writeRaster(bfal_ESMmean_hist, filename = "BFAL_ESMmean_2070-2100.tif", format = "GTiff", overwrite = TRUE)
writeRaster(bfal_ESMsd_hist, filename = "BFAL_ESMsd_2070-2100.tif", format = "GTiff", overwrite = TRUE)

# ∆HSI ####
setwd(ipsl_dir)
hist_ipsl <- raster("Climatol_1985-2015/BFAL_OccProbClimatol_IPSL.tif")
fut_ipsl <- raster("Climatol_2070-2100/BFAL_OccProbClimatol_IPSL.tif")
dif_ipsl_bfal <- fut_ipsl - hist_ipsl
writeRaster(dif_ipsl_bfal, filename = "BFAL_∆HSI_IPSL.tif", format = "GTiff", overwrite = TRUE)

# _____________________ ####
# SEASONAL CLIMATOLOGIES
# bfal ipsl historical (1985-2015) ####
# Example: createa seasonal climatolgies for Species == black-footed albatross
# Earth-system model == IPSL
# time period == Historical (1985-2015)
predDir="bfal/dailyP_IPSL/"  # change to the folder where your daily projections (1980-2100) are stored
setwd(predDir)

fpath <- getwd()
raster_files <- list.files(pattern = "\\.gri$") # get files

extract_date <- function(filename) {
  date_string <- gsub("bfal43_", "", filename) # Remove 'bfal_' prefix
  date_string <- gsub("_mean\\.gri", "", date_string) # Remove '_mean.gri' suffix
  return(ymd(date_string)) # Parse date
}

# Filter files from 1985 to 2015
raster_files <- raster_files[extract_date(raster_files) >= ymd("1985-01-01") & 
                               extract_date(raster_files) <= ymd("2015-12-31")]

seasonal_means_folder <- file.path(fpath, "seasonal_means")
if (!dir.exists(seasonal_means_folder)) {
  dir.create(seasonal_means_folder)
}

seasons <- c(
  "Davidson",
  "Oceanic", "Upwelling")

season_dates <- list(
  Davidson = c("11-15", "03-14"),
  Oceanic = c("08-15", "11-14"),
  Upwelling = c("03-15", "08-15")
)

for (season in seasons) {
  start_year <- 1985
  end_year <- 2015
  
  start_date <- ymd(paste0(start_year, "-", season_dates[[season]][1]))
  end_date <- ymd(paste0(end_year, "-", season_dates[[season]][2]))
  
  # Filter raster files for the current season and years
  season_raster_files <- raster_files[extract_date(raster_files) >= start_date & 
                                        extract_date(raster_files) <= end_date]
  
  start_month <- month(start_date)
  start_day <- day(start_date)
  end_month <- month(end_date)
  end_day <- day(end_date)
  
  # Filter raster files within the season and specific month and day
  if (season == "Davidson") {
    season_raster_files <- season_raster_files[(month(extract_date(season_raster_files)) %in% c(1,2, 12)) |
                                                 (month(extract_date(season_raster_files)) == start_month & 
                                                    day(extract_date(season_raster_files)) >= start_day) |
                                                 (month(extract_date(season_raster_files)) == end_month & 
                                                    day(extract_date(season_raster_files)) <= end_day)
    ]
  } else {
    season_raster_files <- season_raster_files[(month(extract_date(season_raster_files)) == start_month & 
                                                  day(extract_date(season_raster_files)) >= start_day) |
                                                 (month(extract_date(season_raster_files)) == end_month & 
                                                    day(extract_date(season_raster_files)) <= end_day) |
                                                 (month(extract_date(season_raster_files)) > start_month & 
                                                    month(extract_date(season_raster_files)) < end_month)]
  }
  
  raster_stack <- raster::stack(season_raster_files)
  mean_raster <- mean(raster_stack)
  writeRaster(mean_raster, file.path(seasonal_means_folder, paste0(season, "_Historicalmean.tif")), overwrite = TRUE)
}



# bfal ipsl future (2070-2100) ####
setwd(predDir) # directory where daily projections are stored

fpath <- getwd()
raster_files <- list.files(pattern = "\\.gri$") # get files

extract_date <- function(filename) {
  date_string <- gsub("bfal43_", "", filename) # Remove 'bfal_' prefix
  date_string <- gsub("_mean\\.gri", "", date_string) # Remove '_mean.gri' suffix
  return(ymd(date_string)) # Parse date
}

# Filter files from 2070 to 2100
raster_files <- raster_files[extract_date(raster_files) >= ymd("2070-01-01") & 
                               extract_date(raster_files) <= ymd("2100-12-31")]

seasonal_means_folder <- file.path(fpath, "seasonal_means")
if (!dir.exists(seasonal_means_folder)) {
  dir.create(seasonal_means_folder)
}

seasons <- c(
  "Davidson",
  "Oceanic", "Upwelling")

season_dates <- list(
  Davidson = c("11-15", "03-14"),
  Oceanic = c("08-15", "11-14"),
  Upwelling = c("03-15", "08-15")
)

for (season in seasons) {
  start_year <- 2070
  end_year <- 2100
  
  start_date <- ymd(paste0(start_year, "-", season_dates[[season]][1]))
  end_date <- ymd(paste0(end_year, "-", season_dates[[season]][2]))
  
  # Filter raster files for the current season and years
  season_raster_files <- raster_files[extract_date(raster_files) >= start_date & 
                                        extract_date(raster_files) <= end_date]
  
  start_month <- month(start_date)
  start_day <- day(start_date)
  end_month <- month(end_date)
  end_day <- day(end_date)
  
  # Filter raster files within the season and specific month and day
  if (season == "Davidson") {
    season_raster_files <- season_raster_files[(month(extract_date(season_raster_files)) %in% c(1,2, 12)) |
                                                 (month(extract_date(season_raster_files)) == start_month & 
                                                    day(extract_date(season_raster_files)) >= start_day) |
                                                 (month(extract_date(season_raster_files)) == end_month & 
                                                    day(extract_date(season_raster_files)) <= end_day)
    ]
  } else {
    season_raster_files <- season_raster_files[(month(extract_date(season_raster_files)) == start_month & 
                                                  day(extract_date(season_raster_files)) >= start_day) |
                                                 (month(extract_date(season_raster_files)) == end_month & 
                                                    day(extract_date(season_raster_files)) <= end_day) |
                                                 (month(extract_date(season_raster_files)) > start_month & 
                                                    month(extract_date(season_raster_files)) < end_month)]
  }
  
  raster_stack <- raster::stack(season_raster_files)
  mean_raster <- mean(raster_stack)
  writeRaster(mean_raster, file.path(seasonal_means_folder, paste0(season, "_Futuremean.tif")), overwrite = TRUE)
}



# _____________________#####
# TIME SERIES BFAL IPSL ####
# Domain ####
setwd(ipsl.dir)  # directory where IPSL projections for BFAL are stored
rm(list=ls())

raster_files <- list.files(pattern = 'bfal_.*\\.grd$')
raster_list <- list()

for(file in raster_files) {
  date_str <- gsub("bfal_", "", gsub("_mean.grd", '', file))
  date <- as.Date(date_str)
  if (date >= as.Date("1980-01-01") && date <= as.Date('2100-12-31')) {
    raster_list[[date_str]] <- raster(file)
  }
}

occ_prob <- sapply(raster_list, function(r) mean(r[], na.rm = TRUE))

occ_prob_df <- data.frame(
  Date = as.Date(names(occ_prob), format = "%Y-%m-%d"),
  Prob = occ_prob
)

occ_prob_df$Species <- 'BFAL'

write.csv(occ_prob_df, file = 'bfal_MeanOccProb_TimeSeries_1980-2100_DOMAIN.csv')  

# ggplot(occ_prob_df, aes(x = Date, y = Prob)) + geom_line() + geom_smooth()

# to 273 km offshore ####
rm(list=ls())
data.dir <- "~/SeabirdModels/Data" # directory where shore distnace raster is stored
setwd(data.dir)
dist <- raster('Distance_from_shore.grd')

setwd(ipsl.dir) # directory where IPSL projections for BFAL are stored
raster_files <- list.files(pattern = 'bfal_.*\\.grd$')
raster_list <- list()

for(file in raster_files) {
  date_str <- gsub("bfal_", "", gsub("_mean.grd", '', file))
  date <- as.Date(date_str)
  if (date >= as.Date("1980-01-01") && date <= as.Date('2100-12-31')) {
    raster_data <- raster(file)
    raster_data[dist[] >= 273000] <- NA # apply dist mask
    raster_list[[date_str]] <- raster_data
  }
}

occ_prob <- sapply(raster_list, function(r) mean(r[], na.rm = TRUE))

occ_prob_df <- data.frame(
  Date = as.Date(names(occ_prob), format = "%Y-%m-%d"),
  Prob = occ_prob
)

occ_prob_df$Species <- 'BFAL'

write.csv(occ_prob_df, file = 'bfal_MeanOccProb_TimeSeries_1980-2100_ShelfSlope.csv')  


# to 200m (on-shelf) ####
rm(list=ls())

setwd("/Volumes/Triple_Bottom_Line/Data/EcoROMS/ROMS_predict/1981-01-01")
z <- raster('z.grd')

setwd(ipsl.dir)  # directory where IPSL projections for BFAL are stored

raster_files <- list.files(pattern = 'bfal_.*\\.grd$')
raster_list <- list()

for(file in raster_files) {
  date_str <- gsub("bfal_", "", gsub("_mean.grd", '', file))
  date <- as.Date(date_str)
  if (date >= as.Date("1980-01-01") && date <= as.Date('2100-12-31')) {
    raster_data <- raster(file)
    raster_data[z[] <= -200] <- NA # apply depth mask
    raster_list[[date_str]] <- raster_data
  }
}

occ_prob <- sapply(raster_list, function(r) mean(r[], na.rm = TRUE))

occ_prob_df <- data.frame(
  Date = as.Date(names(occ_prob), format = "%Y-%m-%d"),
  Prob = occ_prob
)

occ_prob_df$Species <- 'BFAL'

write.csv(occ_prob_df, file = 'bfal_MeanOccProb_TimeSeries_1980-2100_200m.csv')  

# CA WEAs ####
rm(list=ls())
wind.dir <- '~/SeabirdModels/Shapefiles/WindAreas' # directory where wind area shapefiles are stored
setwd(wind.dir)
list.files()

shape <- read_sf('CA_FinalLeaseAreas_Shapefiles/LeaseAreas_CA_Outlines_2022_09_08.shp')
shape <- st_transform(shape, crs = 4326)
shape <- st_zm(shape, drop = TRUE)
shape_humboldt <- subset(shape, Lease_Numb == 'OCS-P 0561' | Lease_Numb == 'OCS-P 0562')
shape_morro <- subset(shape, Lease_Numb == 'OCS-P 0563' | Lease_Numb == 'OCS-P 0564' | Lease_Numb == 'OCS-P 0565')
shape_humboldt <- st_union(shape_humboldt)
shape_morro <- st_union(shape_morro)

shape <- as(shape, "Spatial") # Convert sf object to SpatialPolygonsDataFrame
shape_humboldt <- as(shape_humboldt, "Spatial") # Convert sf object to SpatialPolygonsDataFrame
shape_morro <- as(shape_morro, "Spatial") # Convert sf object to SpatialPolygonsDataFrame


setwd(ipsl.dir)

raster_files <- list.files(pattern = 'bfal_.*\\.grd$')
raster_list <- list()

for(file in raster_files) {
  date_str <- gsub("bfal_", "", gsub("_mean.grd", '', file))
  date <- as.Date(date_str)
  if (date >= as.Date("1980-01-01") && date <= as.Date('2100-12-31')) {
    raster_data <- raster(file)
    raster_data_cropped <- as.data.frame(raster::extract(raster_data, shape, fun = mean, na.rm = TRUE))
    raster_list[[date_str]] <- raster_data_cropped
  }
}

occ_prob_df <- data.frame(
  Date = as.Date(names(raster_list), format = "%Y-%m-%d"),
  Prob = unlist(raster_list)
)

raster_list <- list()

for(file in raster_files) {
  date_str <- gsub("bfal_", "", gsub("_mean.grd", '', file))
  date <- as.Date(date_str)
  if (date >= as.Date("1980-01-01") && date <= as.Date('2100-12-31')) {
    raster_data <- raster(file)
    raster_data_cropped <- as.data.frame(raster::extract(raster_data, shape_humboldt, fun = mean, na.rm = TRUE))
    raster_list[[date_str]] <- raster_data_cropped
  }
}

occ_prob_df_humboldt <- data.frame(
  Date = as.Date(names(raster_list), format = "%Y-%m-%d"),
  Prob = unlist(raster_list)
)

raster_list <- list()

for(file in raster_files) {
  date_str <- gsub("bfal_", "", gsub("_mean.grd", '', file))
  date <- as.Date(date_str)
  if (date >= as.Date("1980-01-01") && date <= as.Date('2100-12-31')) {
    raster_data <- raster(file)
    raster_data_cropped <- as.data.frame(raster::extract(raster_data, shape_morro, fun = mean, na.rm = TRUE))
    raster_list[[date_str]] <- raster_data_cropped
  }
}

occ_prob_df_morro <- data.frame(
  Date = as.Date(names(raster_list), format = "%Y-%m-%d"),
  Prob = unlist(raster_list)
)

occ_prob_df$area <- row.names(occ_prob_df)
occ_prob_df$lease <- factor(substr(occ_prob_df$area, nchar(occ_prob_df$area), nchar(occ_prob_df$area)))
occ_prob_df$lease2 <-  occ_prob_df$lease

levels(occ_prob_df$lease2) <- unique(shape$Lease_Numb)
levels(occ_prob_df$lease)
levels(occ_prob_df$lease2)

occ_prob_df$Species <- 'BFAL'
occ_prob_df_humboldt$Species <- 'BFAL'
occ_prob_df_morro$Species <- 'BFAL'

# write.csv(occ_prob_df, file = 'bfal_MeanOccProb_TimeSeries_1980-2100_CA_Leases.csv')  
write.csv(occ_prob_df_humboldt, file = 'bfal_MeanOccProb_TimeSeries_1980-2100_HumboldtWEA.csv')  
write.csv(occ_prob_df_morro, file = 'bfal_MeanOccProb_TimeSeries_1980-2100_MorroBayWEA.csv')  

# OR WEAs ####
rm(list=ls())
wind.dir <- '~/SeabirdModels/Shapefiles/WindAreas' # directory where wind area shapefiles are stored
setwd(wind.dir)
list.files()

shape <- read_sf('OR_WEA_GIS_Data/CDR_CDR_OREGON_WEA_2024_BOEM.shp')
shape <- st_transform(shape, crs = 4326)
shape <- st_zm(shape, drop = TRUE)
shape_brook <- subset(shape, Name == 'Brookings')
shape_coos <- subset(shape, Name == 'Coos Bay')

shape_brook <- as(shape_brook, "Spatial") # Convert sf object to SpatialPolygonsDataFrame
shape_coos <- as(shape_coos, "Spatial") # Convert sf object to SpatialPolygonsDataFrame

setwd(ipsl.dir)

raster_files <- list.files(pattern = 'bfal_.*\\.grd$')
raster_list <- list()

for(file in raster_files) {
  date_str <- gsub("bfal_", "", gsub("_mean.grd", '', file))
  date <- as.Date(date_str)
  if (date >= as.Date("1980-01-01") && date <= as.Date('2100-12-31')) {
    raster_data <- raster(file)
    raster_data_cropped <- as.data.frame(raster::extract(raster_data, shape_brook, fun = mean, na.rm = TRUE))
    raster_list[[date_str]] <- raster_data_cropped
  }
}

occ_prob_df <- data.frame(
  Date = as.Date(names(raster_list), format = "%Y-%m-%d"),
  Prob = unlist(raster_list)
)

occ_prob_df$Species <- 'BFAL'
head(occ_prob_df)

write.csv(occ_prob_df, file = 'bfal_MeanOccProb_TimeSeries_1980-2100_BrookingsWEA.csv')  

rm(occ_prob_df)

for(file in raster_files) {
  date_str <- gsub("bfal_", "", gsub("_mean.grd", '', file))
  date <- as.Date(date_str)
  if (date >= as.Date("1980-01-01") && date <= as.Date('2100-12-31')) {
    raster_data <- raster(file)
    raster_data_cropped <- as.data.frame(raster::extract(raster_data, shape_coos, fun = mean, na.rm = TRUE))
    raster_list[[date_str]] <- raster_data_cropped
  }
}

occ_prob_df <- data.frame(
  Date = as.Date(names(raster_list), format = "%Y-%m-%d"),
  Prob = unlist(raster_list)
)


occ_prob_df$Species <- 'BFAL'
head(occ_prob_df)

write.csv(occ_prob_df, file = 'bfal_MeanOccProb_TimeSeries_1980-2100_CoosBayWEA.csv')  

# Monterey Bay NMS ####
rm(list=ls())
shape.dir <- '~/SeabirdModels/Shapefiles' # directory where National Marine Sanctuary (NMS) shapefiles are stored
setwd(shape.dir)
list.files()

shape <- read_sf('Monterey_Bay')
shape <- st_transform(shape, crs = 4326)
shape <- st_zm(shape, drop = TRUE)
shape <- as(shape, "Spatial") # Convert sf object to SpatialPolygonsDataFrame

setwd(ipsl.dir)

raster_files <- list.files(pattern = 'bfal_.*\\.grd$')
raster_list <- list()

for(file in raster_files) {
  date_str <- gsub("bfal_", "", gsub("_mean.grd", '', file))
  date <- as.Date(date_str)
  if (date >= as.Date("1980-01-01") && date <= as.Date('2100-12-31')) {
    raster_data <- raster(file)
    raster_data_cropped <- as.data.frame(raster::extract(raster_data, shape, fun = mean, na.rm = TRUE))
    raster_list[[date_str]] <- raster_data_cropped
  }
}

occ_prob_df <- data.frame(
  Date = as.Date(names(raster_list), format = "%Y-%m-%d"),
  Prob = unlist(raster_list)
)

occ_prob_df$Species <- 'BFAL'
head(occ_prob_df)

write.csv(occ_prob_df, file = 'bfal_MeanOccProb_TimeSeries_1980-2100_MBNMS.csv')  

# Cordell Bank NMS ####
rm(list=ls())
shape.dir <- '~/SeabirdModels/Shapefiles' # directory where National Marine Sanctuary (NMS) shapefiles are stored
setwd(shape.dir)
list.files()

shape <- read_sf('Cordell_Bank')
shape <- st_transform(shape, crs = 4326)
shape <- st_zm(shape, drop = TRUE)
shape <- as(shape, "Spatial") # Convert sf object to SpatialPolygonsDataFrame

setwd(ipsl.dir)

raster_files <- list.files(pattern = 'bfal_.*\\.grd$')
raster_list <- list()

for(file in raster_files) {
  date_str <- gsub("bfal_", "", gsub("_mean.grd", '', file))
  date <- as.Date(date_str)
  if (date >= as.Date("1980-01-01") && date <= as.Date('2100-12-31')) {
    raster_data <- raster(file)
    raster_data_cropped <- as.data.frame(raster::extract(raster_data, shape, fun = mean, na.rm = TRUE))
    raster_list[[date_str]] <- raster_data_cropped
  }
}

occ_prob_df <- data.frame(
  Date = as.Date(names(raster_list), format = "%Y-%m-%d"),
  Prob = unlist(raster_list)
)

occ_prob_df$Species <- 'BFAL'
head(occ_prob_df)

write.csv(occ_prob_df, file = 'bfal_MeanOccProb_TimeSeries_1980-2100_CBNMS.csv')  

# Farallones NMS ####
rm(list=ls())
shape.dir <- '~/SeabirdModels/Shapefiles' # directory where National Marine Sanctuary (NMS) shapefiles are stored
setwd(shape.dir)
list.files()

shape <- read_sf('Greater_Farallones')
shape <- st_transform(shape, crs = 4326)
shape <- st_zm(shape, drop = TRUE)
shape <- as(shape, "Spatial") # Convert sf object to SpatialPolygonsDataFrame

setwd(ipsl.dir)

raster_files <- list.files(pattern = 'bfal_.*\\.grd$')
raster_list <- list()

for(file in raster_files) {
  date_str <- gsub("bfal_", "", gsub("_mean.grd", '', file))
  date <- as.Date(date_str)
  if (date >= as.Date("1980-01-01") && date <= as.Date('2100-12-31')) {
    raster_data <- raster(file)
    raster_data_cropped <- as.data.frame(raster::extract(raster_data, shape, fun = mean, na.rm = TRUE))
    raster_list[[date_str]] <- raster_data_cropped
  }
}

occ_prob_df <- data.frame(
  Date = as.Date(names(raster_list), format = "%Y-%m-%d"),
  Prob = unlist(raster_list)
)

occ_prob_df$Species <- 'BFAL'
head(occ_prob_df)

write.csv(occ_prob_df, file = 'bfal_MeanOccProb_TimeSeries_1980-2100_GFNMS.csv')  

# Chumash Heritage NMS ####
rm(list=ls())
shape.dir <- '~/SeabirdModels/Shapefiles' # directory where National Marine Sanctuary (NMS) shapefiles are stored
setwd(shape.dir)
list.files()

shape <- read_sf('Chumash_Heritage')
shape <- st_transform(shape, crs = 4326)
shape <- st_zm(shape, drop = TRUE)
shape <- as(shape, "Spatial") # Convert sf object to SpatialPolygonsDataFrame

setwd(ipsl.dir)

raster_files <- list.files(pattern = 'bfal_.*\\.grd$')
raster_list <- list()

for(file in raster_files) {
  date_str <- gsub("bfal_", "", gsub("_mean.grd", '', file))
  date <- as.Date(date_str)
  if (date >= as.Date("1980-01-01") && date <= as.Date('2100-12-31')) {
    raster_data <- raster(file)
    raster_data_cropped <- as.data.frame(raster::extract(raster_data, shape, fun = mean, na.rm = TRUE))
    raster_list[[date_str]] <- raster_data_cropped
  }
}

occ_prob_df <- data.frame(
  Date = as.Date(names(raster_list), format = "%Y-%m-%d"),
  Prob = unlist(raster_list)
)

occ_prob_df$Species <- 'BFAL'
head(occ_prob_df)

write.csv(occ_prob_df, file = 'bfal_MeanOccProb_TimeSeries_1980-2100_CHNMS.csv')  

# _____________________ ####
