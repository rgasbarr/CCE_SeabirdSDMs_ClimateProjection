# HEADER ####
# Author: Ryan Gasbarro; rygasbar@ucsc.edu
# Script Purpose: Create maps showing where projected core habitat exists 
# in Wind Energy Areas (WEA) in each time period
# For CCE seabird species using habitat suitability model outputs
# See Gasbarro et al. (2025) PLoS Clim. for information


# Load Climatology raters ####
# *** NOTE *** These rasters were created using the code from script: 5_∆HSI_ClimatolPlots.R
# This Script MUST BE RUN PRIOR TO THIS ONE **************************************

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




# Load WEA data & calculate core habitat ####
wind.dir <- '~/SeabirdModels/Shapefiles/WindAreas' # directory where WEA shapefiles are stored
setwd(wind.dir)

MB <- read_sf('MorroBay_MERGE.shp')
H <- read_sf('Humboldt_MERGE.shp')
B <- read_sf('Brookings_MERGE.shp')
CB <- read_sf('CoosBay_MERGE.shp')

MB <- st_transform(MB, crs = st_crs(sosh_ensem_mean)) 
H  <- st_transform(H,  crs = st_crs(sosh_ensem_mean))
B  <- st_transform(B,  crs = st_crs(sosh_ensem_mean))
CB <- st_transform(CB, crs = st_crs(sosh_ensem_mean))

MB$id <- 'Morro Bay'
H$id <- 'Humboldt'
B$id <- 'Brookings'
CB$id <- 'Coos Bay'

# 2. Combine all WEA polygons
WEAs <- rbind(MB, H, B, CB)
# Re-drop Z after combining polygons
WEAs <- st_zm(WEAs)

# Reproject to match raster (just in case)
WEAs <- st_transform(WEAs, crs = st_crs(sosh_ensem_mean))

# Now convert to Spatial
WEAs_sp <- as(WEAs, "Spatial")
utm_crs <- 26910

# Transform to projected CRS, buffer, then transform back to WGS84
WEAs_proj <- st_transform(WEAs, crs = utm_crs)

# Apply 10 km buffer
WEAs_buffered_proj <- st_buffer(WEAs_proj, dist = 10000)

# Transform back to WGS84 for mapping/raster work
WEAs_buffered <- st_transform(WEAs_buffered_proj, crs = st_crs(sosh_ensem_mean))
WEAs_buffered$id <- WEAs$id

# Get bounding box for each WEA to set individual panel extents
WEA_extents <- st_bbox(WEAs)  # Bounding box for all WEAs

# Create a list to store the raster data for each WEA (with extent) and the corresponding bounding box
raster_list <- list()

sosh_masked <- mask(sosh_ensem_mean, as(WEAs_buffered, "Spatial"))

r_df <- as.data.frame(sosh_masked, xy = TRUE, na.rm = TRUE)
colnames(r_df)[3] <- "sosh"

# Apply the mask to bfal_ensem_mean
bfal_masked <- mask(bfal_ensem_mean, mask_raster, maskvalue = 1)
sosh_masked <- mask(sosh_ensem_mean, mask_raster, maskvalue = 1)
comu_masked <- mask(comu_ensem_mean, mask_raster, maskvalue = 1)
caau_masked <- mask(caau_ensem_mean, mask_raster, maskvalue = 1)
rhau_masked <- mask(rhau_ensem_mean, mask_raster, maskvalue = 1)

bfal_q75 <- quantile(bfal_masked, probs = 0.75, na.rm = TRUE)
comu_q75 <- quantile(comu_masked, probs = 0.75, na.rm = TRUE)
caau_q75 <- quantile(caau_masked, probs = 0.75, na.rm = TRUE)
sosh_q75 <- quantile(sosh_masked, probs = 0.75, na.rm = TRUE)
rhau_q75 <- quantile(rhau_masked, probs = 0.75, na.rm = TRUE)

#Historical
# Loop through each WEA and mask the rasters for all species
for (wea in unique(WEAs$id)) {
  # Mask the rasters for each species by the current WEA
  bfal_masked_wea <- mask(bfal_masked, WEAs_buffered[WEAs_buffered$id == wea, ])
  sosh_masked_wea <- mask(sosh_masked, WEAs_buffered[WEAs_buffered$id == wea, ])
  comu_masked_wea <- mask(comu_masked, WEAs_buffered[WEAs_buffered$id == wea, ])
  caau_masked_wea <- mask(caau_masked, WEAs_buffered[WEAs_buffered$id == wea, ])
  rhau_masked_wea <- mask(rhau_masked, WEAs_buffered[WEAs_buffered$id == wea, ])
  
  # Convert each masked raster to a data.frame
  r_df_wea_bfal <- as.data.frame(bfal_masked_wea, xy = TRUE, na.rm = TRUE)
  r_df_wea_sosh <- as.data.frame(sosh_masked_wea, xy = TRUE, na.rm = TRUE)
  r_df_wea_comu <- as.data.frame(comu_masked_wea, xy = TRUE, na.rm = TRUE)
  r_df_wea_caau <- as.data.frame(caau_masked_wea, xy = TRUE, na.rm = TRUE)
  r_df_wea_rhau <- as.data.frame(rhau_masked_wea, xy = TRUE, na.rm = TRUE)
  
  # Rename columns for each species
  colnames(r_df_wea_bfal)[3] <- "bfal"
  colnames(r_df_wea_sosh)[3] <- "sosh"
  colnames(r_df_wea_comu)[3] <- "comu"
  colnames(r_df_wea_caau)[3] <- "caau"
  colnames(r_df_wea_rhau)[3] <- "rhau"
  
  # Add the 'id' column to each data.frame
  r_df_wea_bfal$id <- wea
  r_df_wea_sosh$id <- wea
  r_df_wea_comu$id <- wea
  r_df_wea_caau$id <- wea
  r_df_wea_rhau$id <- wea
  
  # Combine the data for all  species into one data frame
  r_df_wea <- merge(r_df_wea_bfal, r_df_wea_sosh, by = c("x", "y", "id"))
  r_df_wea <- merge(r_df_wea, r_df_wea_comu, by = c("x", "y", "id"))
  r_df_wea <- merge(r_df_wea, r_df_wea_caau, by = c("x", "y", "id"))
  r_df_wea <- merge(r_df_wea, r_df_wea_rhau, by = c("x", "y", "id"))
  
  # Get the bounding box for the current WEA to set individual limits for each panel
  wea_bbox <- st_bbox(WEAs[WEAs$id == wea, ])
  
  # Store raster data and bounding box together in the list
  raster_list[[wea]] <- list(data = r_df_wea, bbox = wea_bbox)
}

# Combine all raster data into one data frame
r_df_all <- do.call(rbind, lapply(raster_list, function(x) x$data))
r_df_all$period <- 'historical'

# Reshape the data for ggplot (long format)
r_df_long <- tidyr::pivot_longer(r_df_all, cols = c("bfal", "sosh", "comu", "caau", "rhau"), 
                                 names_to = "species", values_to = "value")

r_df_long$species <- factor(r_df_long$species, levels = c("bfal", "sosh", "comu", "caau", "rhau"))
r_df_long$species <- toupper(as.character(r_df_long$species))

# Reorder the factor levels in uppercase
r_df_long$species <- factor(r_df_long$species, levels = c("BFAL", "SOSH", "COMU", "CAAU", "RHAU"))

r_df_long <- r_df_long %>%
  mutate(
    is_core = case_when(
      species == "BFAL" & value > bfal_q75 ~ TRUE,
      species == "SOSH" & value > sosh_q75 ~ TRUE,
      species == "COMU" & value > comu_q75 ~ TRUE,
      species == "CAAU" & value > caau_q75 ~ TRUE,
      species == "RHAU" & value > rhau_q75 ~ TRUE,
      TRUE ~ FALSE
    )
  )

r_df_long$is_core <- as.factor(ifelse(r_df_long$is_core == TRUE, 1, 0))


# Future 
for (wea in unique(WEAs$id)) {
  # Mask the rasters for each species by the current WEA
  bfal_masked_wea <- mask(bfal_masked_fut, WEAs_buffered[WEAs_buffered$id == wea, ])
  sosh_masked_wea <- mask(sosh_masked_fut, WEAs_buffered[WEAs_buffered$id == wea, ])
  comu_masked_wea <- mask(comu_masked_fut, WEAs_buffered[WEAs_buffered$id == wea, ])
  caau_masked_wea <- mask(caau_masked_fut, WEAs_buffered[WEAs_buffered$id == wea, ])
  rhau_masked_wea <- mask(rhau_masked_fut, WEAs_buffered[WEAs_buffered$id == wea, ])
  
  # Convert each masked raster to a data.frame
  r_df_wea_bfal <- as.data.frame(bfal_masked_wea, xy = TRUE, na.rm = TRUE)
  r_df_wea_sosh <- as.data.frame(sosh_masked_wea, xy = TRUE, na.rm = TRUE)
  r_df_wea_comu <- as.data.frame(comu_masked_wea, xy = TRUE, na.rm = TRUE)
  r_df_wea_caau <- as.data.frame(caau_masked_wea, xy = TRUE, na.rm = TRUE)
  r_df_wea_rhau <- as.data.frame(rhau_masked_wea, xy = TRUE, na.rm = TRUE)
  
  # Rename columns for each species
  colnames(r_df_wea_bfal)[3] <- "bfal"
  colnames(r_df_wea_sosh)[3] <- "sosh"
  colnames(r_df_wea_comu)[3] <- "comu"
  colnames(r_df_wea_caau)[3] <- "caau"
  colnames(r_df_wea_rhau)[3] <- "rhau"
  
  # Add the 'id' column to each data.frame
  r_df_wea_bfal$id <- wea
  r_df_wea_sosh$id <- wea
  r_df_wea_comu$id <- wea
  r_df_wea_caau$id <- wea
  r_df_wea_rhau$id <- wea
  
  # Combine the data for all  species into one data frame
  r_df_wea <- merge(r_df_wea_bfal, r_df_wea_sosh, by = c("x", "y", "id"))
  r_df_wea <- merge(r_df_wea, r_df_wea_comu, by = c("x", "y", "id"))
  r_df_wea <- merge(r_df_wea, r_df_wea_caau, by = c("x", "y", "id"))
  r_df_wea <- merge(r_df_wea, r_df_wea_rhau, by = c("x", "y", "id"))
  
  # Get the bounding box for the current WEA to set individual limits for each panel
  wea_bbox <- st_bbox(WEAs[WEAs$id == wea, ])
  
  # Store raster data and bounding box together in the list
  raster_list[[wea]] <- list(data = r_df_wea, bbox = wea_bbox)
}

# Combine all raster data into one data frame
r_df_all_fut <- do.call(rbind, lapply(raster_list, function(x) x$data))
r_df_all_fut$period <- 'future'

# Reshape the data for ggplot (long format)
r_df_long_fut <- tidyr::pivot_longer(r_df_all_fut, cols = c("bfal", "sosh", "comu", "caau", "rhau"), 
                                 names_to = "species", values_to = "value_fut")

r_df_long_fut$species <- factor(r_df_long_fut$species, levels = c("bfal", "sosh", "comu", "caau", "rhau"))
r_df_long_fut$species <- toupper(as.character(r_df_long_fut$species))

# Reorder the factor levels in uppercase
r_df_long_fut$species <- factor(r_df_long_fut$species, levels = c("BFAL", "SOSH", "COMU", "CAAU", "RHAU"))

r_df_long_fut <- r_df_long_fut %>%
  mutate(
    is_core = case_when(
      species == "BFAL" & value_fut > bfal_q75 ~ TRUE,
      species == "SOSH" & value_fut > sosh_q75 ~ TRUE,
      species == "COMU" & value_fut > comu_q75 ~ TRUE,
      species == "CAAU" & value_fut > caau_q75 ~ TRUE,
      species == "RHAU" & value_fut > rhau_q75 ~ TRUE,
      TRUE ~ FALSE
    )
  )

r_df_long_fut$is_core <- as.factor(ifelse(r_df_long_fut$is_core == TRUE, 1, 0))

head(r_df_long)
head(r_df_long_fut)

r_df_long$value_fut <- r_df_long_fut$value_fut
r_df_long$is_core_fut <- r_df_long_fut$is_core

r_df_mb <- subset(r_df_long, id == 'Morro Bay')
r_df_h <- subset(r_df_long, id == 'Humboldt')
r_df_b <- subset(r_df_long, id == 'Brookings')
r_df_cb <- subset(r_df_long, id == 'Coos Bay')

r_df_mb_core <- subset(r_df_mb, is_core == 1)
r_df_h_core <- subset(r_df_h, is_core == 1)
r_df_b_core <- subset(r_df_b, is_core == 1)
r_df_cb_core <- subset(r_df_cb, is_core == 1)

r_df_mb_core_fut <- subset(r_df_mb, is_core_fut == 1)
r_df_h_core_fut <- subset(r_df_h, is_core_fut == 1)
r_df_b_core_fut <- subset(r_df_b, is_core_fut == 1)
r_df_cb_core_fut <- subset(r_df_cb, is_core_fut == 1)

r_df_bfal <- subset(r_df_long, species == 'BFAL')
r_df_sosh <- subset(r_df_long, species == 'SOSH')
r_df_comu <- subset(r_df_long, species == 'COMU')
r_df_caau <- subset(r_df_long, species == 'CAAU')
r_df_rhau <- subset(r_df_long, species == 'RHAU')

mean_labels <- r_df_long %>%
  group_by(id, species) %>%
  summarize(mean_value = mean(value, na.rm = TRUE), .groups = "drop")

# PLOTS ####
a_wea <- ggplot() +
  geom_tile(data = r_df_mb, aes(x = x, y = y, fill = value)) +
  geom_point(data = r_df_mb_core, aes(x = x, y = y), shape = 16, color = "black", size =1.25, alpha = 0.8) +
  geom_point(data = r_df_mb_core_fut, aes(x = x, y = y), shape = 16, color = "white", size =1.25, alpha = 1) +
  geom_sf(data = MB, fill = 'transparent', col = 'black') + 
  facet_wrap(~ species, ncol = 5) +
  scale_fill_viridis_c(name = "HSI", option = "plasma", limits = c(0, 0.75)) +
  theme_minimal() +
  guides(fill = 'none') +
  theme(
    legend.position = 'top',
    strip.text = element_text(size = 10),  # Adjust facet label size if necessary
    axis.text = element_blank()     # Adjust axis label size
  ) +
  labs(
    x = NULL, 
    y = NULL
  )

b_wea <- ggplot() +
  geom_tile(data = r_df_h, aes(x = x, y = y, fill = value)) +
  geom_point(data = r_df_h_core, aes(x = x, y = y), shape = 16, color = "black", size =1.25, alpha = 0.8) +
  geom_point(data = r_df_h_core_fut, aes(x = x, y = y), shape = 16, color = "white", size =1.25, alpha = 1) +
  geom_sf(data = H, fill = 'transparent', col = 'black') + 
  facet_wrap(~ species, ncol = 5) +
  guides(fill = 'none') +
  scale_fill_viridis_c(name = "HSI", option = "plasma", limits = c(0, 0.75)) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10),  # Adjust facet label size if necessary
    axis.text = element_blank()     # Adjust axis label size
  ) +
  labs(
    x = NULL,
    y = NULL
  )

c_wea <- ggplot() +
  geom_tile(data = r_df_b, aes(x = x, y = y, fill = value)) +
  geom_point(data = r_df_b_core, aes(x = x, y = y), shape = 16, color = "black", size =1.25, alpha = 0.8) +
  geom_point(data = r_df_b_core_fut, aes(x = x, y = y), shape = 16, color = "white", size =1.25, alpha = 1) +
  geom_sf(data = B, fill = 'transparent', col = 'black') + 
  facet_wrap(~ species, ncol = 5) +
  scale_fill_viridis_c(name = "HSI", option = "plasma", limits = c(0, 0.75)) +
  theme_minimal() +
  guides(fill = 'none') +
  theme(
    strip.text = element_text(size = 10),  # Adjust facet label size if necessary
    axis.text = element_blank()     # Adjust axis label size
  ) +
  labs(
    x = NULL,
    y = NULL
  )

d_wea <- ggplot() +
  geom_tile(data = r_df_cb, aes(x = x, y = y, fill = value)) +
  geom_point(data = r_df_cb_core, aes(x = x, y = y), shape = 16, color = "black", size =1.25, alpha = 0.8) +
  geom_point(data = r_df_cb_core_fut, aes(x = x, y = y), shape = 16, color = "white", size =1.25, alpha = 1) +
  geom_sf(data = CB, fill = 'transparent', col = 'black') + 
  facet_wrap(~ species, ncol = 5) +
  scale_fill_viridis_c(name = "HSI", option = "plasma", limits = c(0, 0.75)) +
  theme_minimal() +
  # guides(fill = 'none') +
  theme(legend.position = 'top',
    strip.text = element_text(size = 10),  # Adjust facet label size if necessary
    axis.text = element_blank()     # Adjust axis label size
  ) +
  labs(
    x = NULL,
    y = NULL
  ) 

grid.arrange(d_wea, c_wea, b_wea, a_wea, ncol = 1)


# Core habitat richness (# of species with core habitat in a given pixel) plots
mb_rich <- r_df_mb %>%
  dplyr::group_by(x, y, id) %>%
  dplyr::summarise(
    richness_hist = sum(as.numeric(as.character(is_core))),       # historical
    richness_fut = sum(as.numeric(as.character(is_core_fut)))     # future
  ) %>%
  dplyr::ungroup()

mb_rich_long <- mb_rich %>%
  tidyr::pivot_longer(
    cols = dplyr::starts_with("richness"),
    names_to = "period",
    values_to = "richness"
  ) %>%
  dplyr::mutate(
    period = dplyr::recode(
      period,
      richness_hist = "Historical",
      richness_fut = "Future"
    )
  )

h_rich <- r_df_h %>%
  dplyr::group_by(x, y, id) %>%
  dplyr::summarise(
    richness_hist = sum(as.numeric(as.character(is_core))),       # historical
    richness_fut = sum(as.numeric(as.character(is_core_fut)))     # future
  ) %>%
  dplyr::ungroup()

h_rich_long <- h_rich %>%
  tidyr::pivot_longer(
    cols = dplyr::starts_with("richness"),
    names_to = "period",
    values_to = "richness"
  ) %>%
  dplyr::mutate(
    period = dplyr::recode(
      period,
      richness_hist = "Historical",
      richness_fut = "Future"
    )
  )

b_rich <- r_df_b %>%
  dplyr::group_by(x, y, id) %>%
  dplyr::summarise(
    richness_hist = sum(as.numeric(as.character(is_core))),       # historical
    richness_fut = sum(as.numeric(as.character(is_core_fut)))     # future
  ) %>%
  dplyr::ungroup()

b_rich_long <- b_rich %>%
  tidyr::pivot_longer(
    cols = dplyr::starts_with("richness"),
    names_to = "period",
    values_to = "richness"
  ) %>%
  dplyr::mutate(
    period = dplyr::recode(
      period,
      richness_hist = "Historical",
      richness_fut = "Future"
    )
  )

cb_rich <- r_df_cb %>%
  dplyr::group_by(x, y, id) %>%
  dplyr::summarise(
    richness_hist = sum(as.numeric(as.character(is_core))),       # historical
    richness_fut = sum(as.numeric(as.character(is_core_fut)))     # future
  ) %>%
  dplyr::ungroup()

cb_rich_long <- cb_rich %>%
  tidyr::pivot_longer(
    cols = dplyr::starts_with("richness"),
    names_to = "period",
    values_to = "richness"
  ) %>%
  dplyr::mutate(
    period = dplyr::recode(
      period,
      richness_hist = "Historical",
      richness_fut = "Future"
    )
  )

richness_wea <- rbind(mb_rich_long, h_rich_long, b_rich_long, cb_rich_long)
richness_wea$period <- factor(richness_wea$period, levels = c("Historical", "Future"))

pals <- wesanderson::wes_palette('Zissou1', type = 'continuous')

mb_rich_plot <- ggplot() +
  geom_tile(data = subset(richness_wea, id == 'Morro Bay'), aes(x = x, y = y, fill = richness)) +
  geom_sf(data = MB, fill = 'transparent', lwd = 0.75, col = 'black') +
  facet_wrap(id ~ period, ncol = 2) +
  scale_fill_gradientn(colors = pals, limits = c(0,5)) +
  theme_minimal() +
  guides() +
  theme_minimal() +
  theme(
    legend.position = 'none',
    strip.text = element_blank(),
    axis.text = element_blank()     # Adjust axis label size
  ) +
  labs(
    x = NULL, 
    y = NULL
  )

h_rich_plot <- ggplot() +
  geom_tile(data = subset(richness_wea, id == 'Humboldt'), aes(x = x, y = y, fill = richness)) +
  # geom_path(data = mbnms_df, aes(x = X, y = Y, group = L1), color = 'black') +
  geom_sf(data = H, fill = 'transparent', lwd = 0.75, col = 'black') +
  facet_wrap(id ~ period, ncol = 2) +
  scale_fill_gradientn(colors = pals, limits = c(0,5)) +
  theme_minimal() +
  guides() +
  theme_minimal() +
  theme(
    legend.position = 'none',
    strip.text = element_blank(), 
    axis.text = element_blank()     # Adjust axis label size
  ) +
  labs(
    x = NULL, 
    y = NULL
  )

b_rich_plot <- ggplot() +
  geom_tile(data = subset(richness_wea, id == 'Brookings'), aes(x = x, y = y, fill = richness)) +
  # geom_path(data = mbnms_df, aes(x = X, y = Y, group = L1), color = 'black') +
  geom_sf(data = B, fill = 'transparent', lwd = 0.75, col = 'black') +
  facet_wrap(id ~ period, ncol = 2) +
  scale_fill_gradientn(colors = pals, , limits = c(0,5)) +
  theme_minimal() +
  guides() +
  theme_minimal() +
  theme(
    legend.position = 'none',
    strip.text = element_blank(),
    axis.text = element_blank()     # Adjust axis label size
  ) +
  labs(
    x = NULL, 
    y = NULL
  )

cb_rich_plot <- ggplot() +
  geom_tile(data = subset(richness_wea, id == 'Coos Bay'), aes(x = x, y = y, fill = richness)) +
  # geom_path(data = mbnms_df, aes(x = X, y = Y, group = L1), color = 'black') +
  geom_sf(data = CB, fill = 'transparent', lwd = 0.75, col = 'black') +
  facet_wrap(id ~ period, ncol = 2) +
  scale_fill_gradientn(colors = pals, limits = c(0,5)) +
  theme_minimal() +
  guides() +
  theme_minimal() +
  theme(
    legend.position = 'none',
    strip.text = element_blank(),
    axis.text = element_blank()     # Adjust axis label size
  ) +
  labs(
    x = NULL, 
    y = NULL
  )

# grid.arrange(cb_rich_plot, b_rich_plot, h_rich_plot, mb_rich_plot, ncol = 1)








