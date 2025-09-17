# CCE_SeabirdSDMs_ClimateProjection

This repository contains the code for the analyses and workflows associated with:

Gasbarro et al. (in revision)

Authors:
Ryan Gasbarro*, David J. Ainley, Kelly S. Andrews, Lisa T. Ballance, Hannah Blondin, Steven Bograd, Stephanie Brodie, Megan Cimino, Thomas Clay, Aspen Ellis, Joseph Evenson, John C. Field, Elliott L. Hazen, Michael Jacox, Jaime Jahncke, Trevor Joyce, Jeffery B. Leirness, Danielle Lipski, Barbara Muhling, Nerea Lezama-Ochoa, Mercedes Pozo-Buil, Adena Schonfeld, Amanda Warlick, Heather Welch, Jen Zamon, Kelly M. Zilliacus, Jarrod A. Santora

*Please contact Ryan Gasbarro (rygasbar@ucsc.edu) for questions about the code and/or analyses.

## Abstract
Climate-induced changes in ocean conditions are altering the distribution of species' suitable habitat in space and time. Therefore, it is important to understand how species habitat use may change across current management boundaries (e.g. marine protected areas) and identify potential future risks that reduce the effectiveness of fixed boundaries or cause negative interactions between wildlife and human ocean-use sectors. Here, we used presence and absence records from a compilation of > 132,000 ship-based and aerial at-sea visual survey transect segments collected from 1980-2017 to fit species distribution models (SDMs) for five abundant and ecologically important seabird species in the California Current Ecosystem, including both resident (common murre, Cassin’s auklet, and rhinoceros auklet) and seasonal migrant (sooty shearwater, black-footed albatross) species with different life-histories. We then projected their daily habitat suitability from 1980 to 2100 using an ensemble of three dynamically downscaled, high-resolution (0.1º) climate projections for the CCE. We then assessed and compared long-term changes in both mean conditions and intra-annual (seasonal) variability within four National Marine Sanctuaries and four proposed areas for offshore wind energy development (Wind Energy Areas; WEAs) in the CCE. Overall, sea surface temperature, bottom depth, daylength, and biogeographic province were the most important variables with contributions to each SDM being species-specific. All five species displayed a negative relationship with increasing temperatures that was most pronounced in the two auklet species. This led to notable declines in habitat suitability index scores throughout the CCE most pronounced south of Point Conception that emerged from historical variability for all species except sooty shearwater. Despite the overall negative trend in habitat suitability over time, we identified extensive species-specific seasonal refugia that highlight potential changes in the intra–annual occurrence of suitable habitat. The observed changes in suitable habitat suggest that our perception of the conservation benefits of marine sanctuaries and the potential interactions between seabirds and new ocean-use development could be significantly different towards the end of the century, and that many impacts may occur by the middle of the century. Thus, it is important to consider the effects of future ocean conditions on species habitat suitability within marine spatial management and planning processes.

## Table of Contents
Code – Contains scripts detailing the SDM analysis workflow, creation of time-series and climatogy maps from SDM outputs, and analyses and visualizations of time-series data. 

Data to run the scripts is contained in the following linked Zenodo repository:
https://doi.org/10.5281/zenodo.17145095

### 1_TransectSppPresMaps.R                    
Purpose: Plot presence/absence data for five seabird species in the California Current

### 2_FitSDMs.R    
Purpose: Fit boosted regression tree models to seabird presence & environmental data, inspect radar charts of variable importance for each spp.

### 3_Project.R                                
Purpose: Script showing how boosted regression tree models are projected using daily ocean model (ROMS) data. This example uses one species: black-footed albatross (BFAL), and one earth-system model (IPSL).

### 4_CreateClimatologies&TimeSeries.R         
Purpose: Script detailing process of creating habitat suitability climatologies & time-series. This example uses one species: black-footed albatross (BFAL), and one earth-system model (IPSL).

### 5_∆HSI_ClimatolPlots.R                     
Purpose: Plot historical ensemble mean and change in HSI between future (2070-2100) & historical (1985-2015) time periods for each seabird species.

### 6_SeasonalClimatolPlots.R                 
Purpose: Calulate and plot seasonal climatologies of occurrence probability from seabird SDMs.

### 7_HistoricalVsFutureTimeSeriesAvgsByArea.R 
Purpose: Create maps showing where projected core habitat exists in four of California's National Marine Sanctuaries (NMS) in each time period. 

### 8_CoreHabitatByTimePeriod.R                
Purpose: Calculate and plot changes in projected core habitat for seabird species by time period (e.g. historical vs. future).

### 9_TimeSeriesDecomposition&TimeofEmergence.R
Purpose: Example script showing time-series decomposition and time-of-emergence calculations using  model outputs

### 10_MapCoreHabitat_NMS.R                    
Purpose: Create maps showing where projected core habitat exists in National Marine Sanctuaries (NMS) in each time period

### 11_MapCoreHabitat_WEAs.R  
Purpose: Create maps showing where projected core habitat exists in proposed Wind Energy Areas (WEA) in each time period

### 12a_RandomSplit_CrossVal.R
Purpose: Example script for performing random k-fold cross-validation on seabird SDMs

### 12b_SpatiotemporalCrossVal.R
Purpose: Example script for performing random spatial (by biogeographic province) and temporal (by year) cross-validation on seabird SDMs

### 13_PlotSpatiotemporalCrossValResults.R
Purpose: Plot spatiotemporal cross-validation (see Scripts 12a & 12b) Results

## Details of Article
Gasbarro et al. (in revision)

## Disclaimer
This repository is a scientific product and is not official communication of the National Oceanic and Atmospheric Administration, or the United States Department of Commerce. All NOAA GitHub project content is provided on an "as is" basis and the user assumes responsibility for its use. Any claims against the Department of Commerce or Department of Commerce bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by DOC or the United States Government.
