# HEADER ####

# Seabird BRTs Model Evaluation via cross-validation
# Author: Ryan Gasbarro; rygasbar@ucsc.edu
# Purpose: Example script for performing random k-fold cross-validation on seabird SDMs
# See Gasbarro et al. (2025) PLoS Climate for more information

# Load Packages ####
library(caret)
library(mlbench)
library(gbm)
library(tidyverse)
library(dismo)

rm(list = ls())
set.seed(42)

# Load data ####
data.dir <- "~/SeabirdModels/Data" # * CHANGE TO YOUR DATA DIRECTORY
setwd(data.dir)
dat <- read.csv("SeabirdsExampleData.csv") # df with seabird obs, metadata (e.g. date, lat, long, observing platform), matched ROMS (e.g. temperature) & static (e.g. bathymetry, distance-to-shore) data

caau_data <- subset(bird_data, Species == 'caau')
bfal_data <- subset(bird_data, Species == 'bfal')
sosh_data <- subset(bird_data, Species == 'sosh')
comu_data <- subset(bird_data, Species == 'comu')
rhau_data <- subset(bird_data, Species == 'rhau')
rm(bird_data)

# Set Params ####
gbm.x=c("z", "z_sd", "wind_stress", "sst", "sst_sd", "ild", "bv", "EKE", "daylength", 'province', 'season')
gbm.y = 'presAbs'
family = "bernoulli"
lr = 0.001
tc = 3
bf=0.6
tolerance = .001 # adjust as needed so at least 1000 trees are fit 
verbose = TRUE

#### _______________________________________ ####
########
#K-folds evaluation
kfolds_eval <- function(dataInput, family, gbm.x, gbm.y, lr=lr, tc,tolerance.method,tolerance,response){
  DataInput <- dataInput
  DataInput$Kset <- kfold(DataInput,10) #randomly allocate k groups
  counter=1
  if(family=="bernoulli"){
    Evaluations_kfold <- as.data.frame(matrix(data=0,nrow=10,ncol=4))
    colnames(Evaluations_kfold) <- c("k","Deviance","AUC","TSS")
  }
  
  if(family=="poisson"){
    Evaluations_kfold <- as.data.frame(matrix(data=0,nrow=10,ncol=5))
    colnames(Evaluations_kfold) <- c("pearson","spearman","Deviance","RMSE","AVE")
  }
  
  for (k in 1:10){
    print(k)
    DataInput_train <- DataInput[DataInput$Kset!=k,]
    DataInput_test <- DataInput[DataInput$Kset==k,]
    DataInput.kfolds <- gbm.step(data=DataInput_train, gbm.x= gbm.x, gbm.y = gbm.y, 
                                 family=family, tree.complexity=tc,
                                 learning.rate = lr, bag.fraction = 0.6, n.trees = 1000, step.size = 50, max.trees = 20000, tolerance.method = tolerance.method,tolerance = tolerance)
    preds <- predict.gbm(DataInput.kfolds, DataInput_test,
                         n.trees=DataInput.kfolds$gbm.call$best.trees, type="response")
    position=grep(response,names(dataInput))
    dev <- calc.deviance(obs=DataInput_test[[response]], pred=preds, calc.mean=TRUE,family = family)
    d <- cbind(DataInput_test[[response]], preds)
    if(family=="bernoulli"){
      pres <- d[d[,1]==1,2]
      abs <- d[d[,1]==0,2]
      e <- evaluate(p=pres, a=abs)
      Evaluations_kfold[counter,1] <- k
      Evaluations_kfold[counter,2] <- dev
      Evaluations_kfold[counter,3] <- e@auc
      Evaluations_kfold[counter,4] <- max(e@TPR + e@TNR-1)
    }
    if(family=="poisson"){
      colnames(d) <- c("observed","predicted")
      d <- as.data.frame(d)
      pear=cor(d$predicted, d$observed, use="na.or.complete",method = "pearson")
      spear=cor(d$predicted, d$observed, use="na.or.complete",method = "spearman")
      # lm(d$observed~d$predicted)
      rmse=sqrt(mean((d$observed - d$predicted)^2)) #RMSE
      ave=mean(d$observed - d$predicted) #AVE
      
      Evaluations_kfold[counter,1] <- pear
      Evaluations_kfold[counter,2] <- spear
      Evaluations_kfold[counter,3] <- dev
      Evaluations_kfold[counter,4] <- rmse
      Evaluations_kfold[counter,5] <- ave
    }
    counter=counter+1 
  }
  return(Evaluations_kfold)}

# DEMO
Species <- "bfal"
bfal_kfold <- kfolds_eval(dataInput = bfal_data, family = family, gbm.x = gbm.x, gbm.y = gbm.y, tolerance = 0.001, tolerance.method = 'auto', lr= lr, tc=3, response = gbm.y)

Species <- "sosh"
sosh_kfold <- kfolds_eval(dataInput = sosh_data, family = family, gbm.x = gbm.x, gbm.y = gbm.y, tolerance = 0.001, tolerance.method = 'auto', lr= lr, tc=3, response = gbm.y)

Species <- "caau"
caau_kfold <- kfolds_eval(dataInput = caau_data, family = family, gbm.x = gbm.x, gbm.y = gbm.y, tolerance = 0.001, tolerance.method = 'auto', lr= lr, tc=3, response = gbm.y)

Species <- "comu"
comu_kfold <- kfolds_eval(dataInput = comu_data, family = family, gbm.x = gbm.x, gbm.y = gbm.y, tolerance = 0.001, tolerance.method = 'auto', lr= lr, tc=3, response = gbm.y)

Species <- "rhau"
rhau_kfold <- kfolds_eval(dataInput = rhau_data, family = family, gbm.x = gbm.x, gbm.y = gbm.y, tolerance = 0.001, tolerance.method = 'auto', lr= lr, tc=3, response = gbm.y)

# Export Cross-val results:

# write.csv(bfal_kfold, file = 'BFAL_10foldCV.csv')
# write.csv(comu_kfold, file = 'COMU_10foldCV.csv')
# write.csv(rhau_kfold, file = 'RHAU_10foldCV.csv')
# write.csv(caau_kfold, file = 'CAAU_10foldCV.csv')
# write.csv(sosh_kfold, file = 'SOSH_10foldCV.csv')


