# HEADER ####
# Example of Seabird SDM Evaluation via leave-one-year-out (LOYO)
# and leave-one-province out (LOPO) cross-validation
# Author: Ryan Gasbarro; rygasbar@ucsc.edu
# Purpose: Example script for performing LOYO cross-validation on seabird SDMs
# See Gasbarro et al. (2025) PLoS Climate for more information


# LOAD PACKAGES ####
library(tidyverse)
library(gridExtra)
library(foreach)
library(dismo)
library(sf)
library(ggsci)
library(ggthemes)
library(adehabitatHR)
library(rgeos)
library(maptools)
# remotes::install_github("eeholmes/WRAP")
library(WRAP)
library(sp)
library(caret)
library(mlbench)
library(lubridate)
library(timeDate)

# Load data ####
data.dir <- "~/SeabirdModels/Data" # * CHANGE TO YOUR DATA DIRECTORY
setwd(data.dir)
data <- read.csv("SeabirdsExampleData.csv") # df with seabird obs, metadata (e.g. date, lat, long, observing platform), matched ROMS (e.g. temperature) & static (e.g. bathymetry, distance-to-shore) data


gbm.x=c("z", "z_sd", "wind_stress", "sst", "sst_sd", "ild", "bv", "EKE", "daylength", 'province', 'season')
gbm.y="presAbs"
response = 'presAbs'
tc = 3
lr = 0.001 # adjust as needed


data$province <- as.factor(data$province)
data$season <- as.factor(data$season)
data$year <- as.numeric(data$year)

data_bfal <- subset(data, Species == 'bfal')
data_caau <- subset(data, Species == 'caau')
data_comu <- subset(data, Species == 'comu')
data_rhau <- subset(data, Species == 'rhau')
data_sosh <- subset(data, Species == 'sosh')

#######
# LOYO CrossVal ####
LOYO_eval <- function(DataInput, gbm.x, gbm.y, lr=lr, tc,tolerance.method,tolerance,family,response){
  DataInput$Year <- as.numeric(DataInput$year)
  Evaluations_LOO <- as.data.frame(matrix(data=0,nrow=1,ncol=4))
  colnames(Evaluations_LOO) <- c("k","Deviance","AUC","TSS")
  counter=1
  for (y in min(DataInput$Year):max(DataInput$Year)){
    print(y)
    DataInput_train <- DataInput[DataInput$Year!=y,]
    DataInput_test <- DataInput[DataInput$Year==y,]
    loyo.gbm <- gbm.step(data=DataInput_train, gbm.x= gbm.x, gbm.y = gbm.y,
                               family='bernoulli', tree.complexity=tc,
                               learning.rate = lr, bag.fraction = 0.6, n.trees = 1000, step.size = 50,tolerance.method = tolerance.method, tolerance = tolerance, max.trees = 20000)
    preds <- predict.gbm(loyo.gbm, DataInput_test,
                         n.trees=loyo.gbm$gbm.call$best.trees, type="response")
    position=grep(response,names(DataInput))
    dev <- calc.deviance(obs=DataInput_test[,position], pred=preds, calc.mean=TRUE)
    d <- cbind(DataInput_test$presAbs, preds)
    pres <- d[d[,1]==1,2]
    abs <- d[d[,1]==0,2]
    if(length(pres)>0 & length(abs)>0){
    e <- evaluate(p=pres, a=abs)

    Evaluations_LOO[counter,1] <- y
    Evaluations_LOO[counter,2] <- dev
    Evaluations_LOO[counter,3] <- e@auc
    Evaluations_LOO[counter,4] <- max(e@TPR + e@TNR-1)
    counter=counter+1
  }
  }
  return(Evaluations_LOO)}

bfal_loyo <- LOYO_eval(data_bfal, gbm.x, gbm.y , lr=0.001, tc=3, tolerance.method = 'auto', tolerance = 0.001, family = 'bernoulli', response = 'presAbs')
# write.csv(bfal_loyo, file = 'BFAL_EvalLOYO.csv')

# caau_loyo <- LOYO_eval(data_caau, gbm.x, gbm.y , lr=0.001, tc=3, tolerance.method = 'auto', tolerance = 0.001, family = 'bernoulli', response = 'presAbs')
write.csv(caau_loyo, file = 'CAAU_EvalLOYO.csv')

comu_loyo <- LOYO_eval(data_comu, gbm.x, gbm.y , lr=0.001, tc=3, tolerance.method = 'auto', tolerance = 0.001, family = 'bernoulli', response = 'presAbs')
# write.csv(comu_loyo, file = 'COMU_EvalLOYO.csv')

rhau_loyo <- LOYO_eval(data_rhau, gbm.x, gbm.y , lr=0.001, tc=3, tolerance.method = 'auto', tolerance = 0.001, family = 'bernoulli', response = 'presAbs')
# write.csv(rhau_loyo, file = 'RHAU_EvalLOYO.csv')

sosh_loyo <- LOYO_eval(data_sosh, gbm.x, gbm.y , lr=0.001, tc=3, tolerance.method = 'auto', tolerance = 0.001, family = 'bernoulli', response = 'presAbs')
# write.csv(sosh_loyo, file = 'SOSH_EvalLOYO.csv')





#######
# LOPO CrossVal ####
LOPO_eval <- function(DataInput, gbm.x, gbm.y, lr=lr, tc, tolerance.method, tolerance, family, response) {
  Evaluations_LOO <- data.frame(k = integer(), Deviance = numeric(), AUC = numeric(), TSS = numeric())
  provinces <- unique(DataInput$province)
  for (province in provinces) {
    print(province)
    # test/train data splits
    DataInput_train <- DataInput[DataInput$province != province, ]
    DataInput_test <- DataInput[DataInput$province == province, ]
    #gbm fitting
    lopo.gbm <- gbm.step(data = DataInput_train, gbm.x = gbm.x, gbm.y = gbm.y,
                         family = 'bernoulli', tree.complexity = tc,
                         learning.rate = lr, bag.fraction = 0.6, n.trees = 8000, 
                         step.size = 50, tolerance.method = tolerance.method, 
                         tolerance = tolerance, max.trees = 20000)
    #make preds
    preds <- predict.gbm(lopo.gbm, DataInput_test, 
                         n.trees = lopo.gbm$gbm.call$best.trees, type = "response")
    position = grep(response, names(DataInput))
    dev <- calc.deviance(obs = DataInput_test[, position], pred = preds, calc.mean = TRUE)
    d <- cbind(DataInput_test$presAbs, preds)
    pres <- d[d[, 1] == 1, 2]
    abs <- d[d[, 1] == 0, 2]
    if (length(pres) > 0 & length(abs) > 0) {
      e <- evaluate(p = pres, a = abs)
      # Append results 
      Evaluations_LOO <- rbind(Evaluations_LOO, data.frame(k = province, Deviance = dev, AUC = e@auc, TSS = max(e@TPR + e@TNR - 1)))
    }
  }
  
  return(Evaluations_LOO)
}


bfal_lopo <- LOPO_eval(data_bfal, gbm.x, gbm.y , lr=0.0025, tc=3, tolerance.method = 'auto', tolerance = 0.001, family = 'bernoulli', response = 'presAbs')
# write.csv(bfal_lopo, file = 'BFAL_Eval_LOPO.csv')
caau_lopo <- LOPO_eval(data_caau, gbm.x, gbm.y , lr=0.0025, tc=3, tolerance.method = 'auto', tolerance = 0.001, family = 'bernoulli', response = 'presAbs')
# write.csv(caau_lopo, file = 'CAAU_Eval_LOPO.csv')
comu_lopo <- LOPO_eval(data_comu, gbm.x, gbm.y , lr=0.0025, tc=3, tolerance.method = 'auto', tolerance = 0.001, family = 'bernoulli', response = 'presAbs')
# write.csv(comu_lopo, file = 'COMU_EvalLOPO.csv')
rhau_lopo <- LOPO_eval(data_rhau, gbm.x, gbm.y , lr=0.0025, tc=3, tolerance.method = 'auto', tolerance = 0.001, family = 'bernoulli', response = 'presAbs')
# write.csv(rhau_lopo, file = 'RHAU_EvalLOPO.csv')
sosh_lopo <- LOPO_eval(data_sosh, gbm.x, gbm.y , lr=0.0025, tc=3, tolerance.method = 'auto', tolerance = 0.001, family = 'bernoulli', response = 'presAbs')
# write.csv(sosh_lopo, file = 'SOSH_EvalLOPO.csv')

