# HEADER ####
# Author: Ryan Gasbarro; rygasbar@ucsc.edu
# Purpose: Example script showing how BRT models are projected using daily ocean model data 
# This example uses one species: black-footed albatross (BFAL)
# And one earth-system model (IPSL)
# See Gasbarro et al. (2025) PLoS Climate for more information

# LOAD PACKAGES ####
rm(list=ls())

library(rlist)
## load libraries
library(dismo)
# library(rgdal)
# install.packages("maptools", repos = "https://packagemanager.posit.co/cran/2023-10-13")
library(maptools)
library(maps)
library(mapdata)
library(raster)
library(ncdf4)
library(dplyr)
library(geosphere)
library(lubridate)
library(terra)
library(sf)
library(rasterVis)

# PROJECT ####
# directory locations ***update as needed to match your organizational structure
droppath="/Volumes/Triple_Bottom_Line/Ryan_working/" 
outDir="/Volumes/Triple_Bottom_Line/Ryan_working/Seabirds_Projections/" 
predDir=paste0(outDir,"bfal/dailyP_IPSL/") 

staticDir=paste0(outDir)
studyarea= readShapeSpatial(paste0(staticDir,"sa_square_coast3.shp"))
template=raster("/Volumes/Triple_Bottom_Line/Ryan_working/Seabirds_Predictions/template.grd")
list.files(paste0(droppath))
list.files(paste0(outDir))

Species="BFAL"
modrep=readRDS("/Volumes/Triple_Bottom_Line/Ryan_working/Seabirds_Projections/p44/mod1_bfal_bernoulli_step_lr0.001_tc3_bf0.6_tol0.001_bernoulli_bfal.rds") # 1 unique models
spname="bfal"
summary(modrep) 

predCIs_ROMS<-function(get_date,spname,modrep,stack,template,outDir,studyarea,droppath){
  
  spDir=paste0(outDir,spname,"/")
  if(!file.exists(spDir)){dir.create(spDir)}
  
  # stack=create_ROMS_daily_stack(get_date = get_date, template=template,staticDir=staticDir, predDir = predDir)
  stack=as.data.frame(stackfile[[i]],stringsAsFactors=F) ## predicts on data.frame, not on a stack
  colnames(stack) <- c("bv", "EKE", "ild", "lat", 'sst_sd', "sst", "sustr", "svstr", "z", 'z_sd', 'wind_stress', 'daylength', 'jday', 'province', 'season')
  desired_order <- c("bv", "EKE", "ild", "lat", 'sst_sd', "sst", "sustr", "svstr", "z", 'z_sd', 'wind_stress', 'daylength', 'jday', 'province', 'season') # ordering to match file structure
  head(stack)
  stack <- stack[, desired_order, drop = FALSE]
  mod_pred <- predict.gbm(modrep, newdata=stack,  type='response')
  
  ## make rasters 
  meanPredR <- setValues(template,mod_pred)%>%mask(.,studyarea)

  ## write rasters 
  writeRaster(meanPredR,paste0(predDir,spname,"_",get_date,"_mean"),overwrite=T)

  
} 

# Daily rasters stored in the following directories:
setwd("/Volumes/Triple_Bottom_Line/Nerea_working/HMS_projections/ROMS_projections/dailyP_IPSL")   
# setwd("/Volumes/Triple_Bottom_Line/Nerea_working/HMS_projections/ROMS_projections/dailyP_GFDL")   
#setwd("/Volumes/Triple_Bottom_Line/Nerea_working//dailyP_HAD/")  
get_date<-list.dirs(getwd(), recursive = F, full.names=F)
get_date<-list.dirs(getwd(), recursive = F, full.names=F)

# function to calculate wind_stress
calculate_wind_stress <- function(sustr, svstr) {
  wind_stress <- sqrt(sustr^2 + svstr^2)
  return(wind_stress)
}

flist <- list()
stackfile <- list()
get_date[1]

stackfile <- list()

# PROJECT ####
for (i in 1:length(get_date)){  
  flist[[i]] <-  list.files(get_date[i], recursive = TRUE, full.names = T, pattern = "z_.1.grd|zsd_.3.grd|svstr.grd|sst.grd|sst_sd.grd|mlat.grd|sustr.grd|ild.grd|bf.grd|EKE.grd$")  
  
  stackfile[[i]] <- stack(flist[[i]])
  names(stackfile[[i]]) <- c('bv', 'EKE', 'ild', 'lat', 'sst_sd', 'sst', 'sustr', 'svstr', 'z', 'z_sd')
  
  # Extract lat from  stack 
  latitude <- as.vector(stackfile[[i]]$lat)
  
  # Calculate Julian day from 'get_date'
  jday <- as.numeric(format(as.Date(get_date[i], format = "%Y-%m-%d"), "%j"))
  
  # assign oceanographic season
  month <- month(get_date[i])
  day <- day(get_date[i])
  
  season <- ifelse((month == 3 & day >= 15) | (month > 3 & month < 8) | (month == 8 & day <= 14), "Upwelling",
                   ifelse((month == 8 & day >= 15) | (month > 8 & month < 11) | (month == 11 & day <= 14), "Oceanic",
                          "Davidson"))
  season <- as.factor(season)
  
  # assign province
  province <- factor()
  province <- ifelse(latitude > 40.4401, "NORTH", NA)
  province <- ifelse(latitude < 40.4401 & latitude > 34.4486, "CENTRAL", province)
  province <- ifelse(latitude < 34.4486, "SOUTH", province)
  province <- as.factor(province)  
  
  # Calculate daylength using geosphere::daylength
  daylength <- geosphere::daylength(latitude, jday)
  daylength_raster <- stackfile[[i]]$lat
  jday_raster <- stackfile[[i]]$lat
  jday_raster[] <- jday
  province_raster <- stackfile[[i]]$lat
  province_raster[] <- province
  names(province_raster) <- 'province'
  season_raster <- stackfile[[i]]$lat
  season_raster[] <- season
  names(season_raster) <- 'season'
  values(daylength_raster) <- daylength
  names(daylength_raster) <- 'daylength'
  names(jday_raster) <- 'jday'
  
  # Calculate wind_stress
  sustr <- as.vector(stackfile[[i]]$sustr)  # Assuming sustr and svstr are named accordingly in the stack
  svstr <- as.vector(stackfile[[i]]$svstr)
  wind_stress <- calculate_wind_stress(sustr, svstr)
  wind_stress_raster <-  stackfile[[i]]$lat
  values(wind_stress_raster) <- wind_stress
  names(wind_stress_raster) <- 'wind_stress'
  
  # Add 'wind_stress' & 'daylength' to the stack
  stackfile[[i]]$wind_stress <- wind_stress_raster
  stackfile[[i]]$daylength <- daylength_raster
  stackfile[[i]]$jday <- jday_raster
  stackfile[[i]]$province <- province_raster
  stackfile[[i]]$season <- season_raster
  
  predCIs_ROMS(get_date = get_date[i],modrep=modrep,spname = spname,studyarea=studyarea,stack=stack,template=template,outDir = predDir,droppath=droppath)
  
  
}

