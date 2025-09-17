# HEADER ####
# Author: Ryan Gasbarro; rygasbar@ucsc.edu
# Purpose: Fit boosted regression tree models to seabird presence & environmental data, inspect radar charts of variable importance for each spp.
# See Gasbarro et al. (2025) PLoS Climate for more information

# LOAD PACKAGES ####
library(tidyverse)
library(gridExtra)
library(foreach)
library(dismo)
# devtools::install_github("ricardo-bion/ggradar")
library(ggradar)
library(sf)
library(ggsci)
library(ggthemes)
library(sp)
# devtools::install_github("JBjouffray/ggBRT")
library(ggBRT)
library(wesanderson)

# LOAD DATA ####
rm(list = ls())
set.seed(42)
data.dir <- "~/SeabirdModels/Data" # * CHANGE TO YOUR DATA DIRECTORY
setwd(data.dir)
dat <- read.csv("SeabirdsExampleData.csv") # df with seabird obs, metadata (e.g. date, lat, long, observing platform), matched ROMS (e.g. temperature) & static (e.g. bathymetry, distance-to-shore) data

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

dat <- subset(dat, lat > 30 & lat <48 & lon > -134 & lon < -115.5) # remove points not on ROMS grid
dat$jday <- lubridate::yday(dat$date)

# Select the columns with spp counts
count_cols1 <- dat[, 15:19]

# Convert count variables to presence-absence format
presence_absence_cols1 <- ifelse(count_cols1 > 0, 1, 0)

dat[, 15:19] <- presence_absence_cols1

dat_pa <- gather(dat, key = "Species", value = "presAbs", 15:19)

dat_pa$season <- ordered(dat_pa$season, levels = c("Davidson", "Upwelling", "Oceanic")) # order the season var


#____________________####
# PRES/ABS BRT MODELS####
set.seed(42)

num_cores <- detectCores()
registerDoParallel(num_cores)

sp = dat_pa$Species %>% unique()

system.time(print(
  foreach(i=1:length(sp),.packages = c("gbm","glue"),.verbose=T,.export = c("dat")) %dopar% {
    print(sp[i])
    
    tryCatch(
      expr ={
        
        main= dat_pa %>% filter(Species==sp[i]) %>% 
          dplyr::select(c(uniqueID, z, z_sd, Species, wind_stress, sst, sst_sd, 
                          ild, bv, ssh, ssh_sd,EKE, daylength, jday, lon, lat, presAbs, province, season)) 
        main=main %>% mutate(presAbs=as.integer(presAbs))
        
        #full model 
        gbm.x=c("z", "z_sd", "wind_stress", "sst", "sst_sd", "ild", "bv", "EKE", "daylength", 'province', 'season')
        family="bernoulli"
        lr=0.001
        tc=3
        bf=0.6
        tolerance = .001
        verbose = TRUE
        type=glue("{sp[i]}")
        
        
        tryCatch(
          expr ={
            brt_step = gbm.step(main,gbm.x=gbm.x,gbm.y="presAbs",family=family,learning.rate = lr, tree.complexity =tc, bag.fraction = bf, max.trees = 20000)
            
            name=glue("model_outs/mod1_{sp[i]}_step_lr{lr}_tc{tc}_bf{bf}_tol{tolerance}_{family}.rds")
            write_rds(brt_step,name)
          },
          error = function(e){
            message(glue("Model not working"))
            print(e)
          }
        )
      },
      error = function(e){
        message(glue("No data"))
        print(e)
      }
    )
  }
))

#____________________####
# EXAMINE MODEL OBJECTS ####
# pres/abs varimport data ####
# rm(list = ls())      # optional to clean up wd
setwd(data.dir)
fold_path <- "model_outs"

dat_pa_bfal <- subset(dat_pa, Species == 'bfal')
dat_pa_rhau <- subset(dat_pa, Species == 'rhau')
dat_pa_caau <- subset(dat_pa, Species == 'caau')
dat_pa_comu <- subset(dat_pa, Species == 'comu')
dat_pa_sosh <- subset(dat_pa, Species == 'sosh')

rds_files <- list.files(fold_path, pattern = "\\.rds$", full.names = TRUE)

# Load each .rds file and store it as an object in the global environment
for (file in rds_files) {
  # Extract the recognizable name from the file name
  object_name <- substr(basename(file), 1, 9)
  
  # Load the .rds file and assign it to the recognizable name
  assign(object_name, readRDS(file))
}


# format model 1 data
model1_df_migrators <- rbind(data.frame(summary(mod1_bfal), sp = rep("bfal", nrow(summary(mod1_bfal)))), 
                   data.frame(summary(mod1_sosh), sp = rep("sosh", nrow(summary(mod1_bfal)))))

model1_df_residents <- rbind(
                   data.frame(summary(mod1_comu), sp = rep("comu", nrow(summary(mod1_bfal)))), 
                   data.frame(summary(mod1_rhau), sp = rep("rhau", nrow(summary(mod1_bfal)))), 
                   data.frame(summary(mod1_caau), sp = rep("caau", nrow(summary(mod1_bfal)))))

model1_df_migrators$rel.inf <- model1_df_migrators$rel.inf/100 #convert to 0-1 scale
model1_df_residents$rel.inf <- model1_df_residents$rel.inf/100 #convert to 0-1 scale

model1_df_migrators <- spread(model1_df_migrators, key = var, value = rel.inf) # convert to wide format for radar plotting
model1_df_residents <- spread(model1_df_residents, key = var, value = rel.inf) # convert to wide format for radar plotting

# radar plots ####
pal1 <- wes_palette("Royal1")
pal2 <- wes_palette("FantasticFox1")

migrators_radar <- ggradar::ggradar(model1_df_migrators, legend.position = 'bottom', legend.title = 'species', 
                 legend.text.size = 12, background.circle.colour = "white", 
                 group.point.size = 2, axis.label.size = 4, grid.label.size = 3) + scale_colour_manual(values = pal1)

residents_radar <- ggradar::ggradar(model1_df_residents, legend.position = 'bottom', legend.title = 'species', 
                                    legend.text.size = 12, background.circle.colour = "white", 
                                    group.point.size = 2, axis.label.size = 4, grid.label.size = 3) + scale_colour_manual(values = pal2)


grid.arrange(migrators_radar, residents_radar, ncol = 2)



