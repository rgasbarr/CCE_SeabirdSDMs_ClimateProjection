# HEADER ####
# Author: Ryan Gasbarro; rygasbar@ucsc.edu
# Script Purpose: Create maps showing where projected core habitat exists 
# in National Marine Sanctuaries (NMS) in each time period
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



# Load NMS data & calculate core habitat ####
shape.dir <- '~/SeabirdModels/Shapefiles' # directory where National Marine Sanctuary (NMS) shapefiles are stored
setwd(shape.dir)

chnms <- read_sf('CHNMS_Boundary_10152024')
mbnms <- read_sf('Monterey_Bay')
cbnms <- read_sf('Cordell_Bank')
gfnms <- read_sf('Greater_Farallones')


chnms <- st_transform(chnms, crs = st_crs(sosh_ensem_mean))
mbnms  <- st_transform(mbnms,  crs = st_crs(sosh_ensem_mean))
cbnms  <- st_transform(cbnms,  crs = st_crs(sosh_ensem_mean))
gfnms <- st_transform(gfnms, crs = st_crs(sosh_ensem_mean))

chnms$id <- 'Chumash Heritage'
mbnms$id <- 'Monterey Bay'
cbnms$id <- 'Cordell Bank'
gfnms$id <- 'Greater Farallones'

chnms <- chnms[, c("id", "geometry")]
mbnms <- mbnms[, c("id", "geometry")]
cbnms <- cbnms[, c("id", "geometry")]
gfnms <- gfnms[, c("id", "geometry")]

# Combine all NMS polygons
NMS <- rbind(chnms, mbnms, cbnms, gfnms)

# Re-drop Z after combining polygons
NMS <- st_zm(NMS)

# Reproject to match raster (just in case)
NMS <- st_transform(NMS, crs = st_crs(sosh_ensem_mean))

# Now convert to Spatial
NMS_sp <- as(NMS, "Spatial")

utm_crs <- 26910

# Transform to projected CRS, buffer, then transform back to WGS84
NMS_proj <- st_transform(NMS, crs = utm_crs)

# Apply 10 km buffer (10,000 meters)
NMS_buffered_proj <- st_buffer(NMS_proj, dist = 10000)

# Transform back to WGS84 for mapping/raster work
NMS_buffered <- st_transform(NMS_buffered_proj, crs = st_crs(sosh_ensem_mean))
NMS_buffered$id <- NMS$id

# Get bounding box for each WEA to set individual panel extents
NMS_extents <- st_bbox(NMS)  # Bounding box for all NMS

# Create a list to store the raster data for each WEA (with extent) and the corresponding bounding box
raster_list <- list()

mask_raster <- dist > 273000

# Apply the mask to bfal_ensem_mean
bfal_masked <- mask(bfal_ensem_mean, mask_raster, maskvalue = 1)
sosh_masked <- mask(sosh_ensem_mean, mask_raster, maskvalue = 1)
comu_masked <- mask(comu_ensem_mean, mask_raster, maskvalue = 1)
caau_masked <- mask(caau_ensem_mean, mask_raster, maskvalue = 1)
rhau_masked <- mask(rhau_ensem_mean, mask_raster, maskvalue = 1)

bfal_masked_fut <- mask(bfal_ensem_mean_fut, mask_raster, maskvalue = 1)
sosh_masked_fut <- mask(sosh_ensem_mean_fut, mask_raster, maskvalue = 1)
comu_masked_fut <- mask(comu_ensem_mean_fut, mask_raster, maskvalue = 1)
caau_masked_fut <- mask(caau_ensem_mean_fut, mask_raster, maskvalue = 1)
rhau_masked_fut <- mask(rhau_ensem_mean_fut, mask_raster, maskvalue = 1)

# Historical
# Loop through each WEA and mask the rasters for all species
for (nms in unique(NMS$id)) {
  # Mask the rasters for each species by the current WEA
  bfal_masked_nms <- mask(bfal_masked, NMS_buffered[NMS_buffered$id == nms, ])
  sosh_masked_nms <- mask(sosh_masked, NMS_buffered[NMS_buffered$id == nms, ])
  comu_masked_nms <- mask(comu_masked, NMS_buffered[NMS_buffered$id == nms, ])
  caau_masked_nms <- mask(caau_masked, NMS_buffered[NMS_buffered$id == nms, ])
  rhau_masked_nms <- mask(rhau_masked, NMS_buffered[NMS_buffered$id == nms, ])
  
  # Convert each masked raster to a data.frame
  r_df_nms_bfal <- as.data.frame(bfal_masked_nms, xy = TRUE, na.rm = TRUE)
  r_df_nms_sosh <- as.data.frame(sosh_masked_nms, xy = TRUE, na.rm = TRUE)
  r_df_nms_comu <- as.data.frame(comu_masked_nms, xy = TRUE, na.rm = TRUE)
  r_df_nms_caau <- as.data.frame(caau_masked_nms, xy = TRUE, na.rm = TRUE)
  r_df_nms_rhau <- as.data.frame(rhau_masked_nms, xy = TRUE, na.rm = TRUE)
  
  # Rename columns for each species
  colnames(r_df_nms_bfal)[3] <- "bfal"
  colnames(r_df_nms_sosh)[3] <- "sosh"
  colnames(r_df_nms_comu)[3] <- "comu"
  colnames(r_df_nms_caau)[3] <- "caau"
  colnames(r_df_nms_rhau)[3] <- "rhau"
  
  # Add the 'id' column to each data.frame
  r_df_nms_bfal$id <- nms
  r_df_nms_sosh$id <- nms
  r_df_nms_comu$id <- nms
  r_df_nms_caau$id <- nms
  r_df_nms_rhau$id <- nms
  
  # Combine the data for all  species into one data frame
  r_df_nms <- merge(r_df_nms_bfal, r_df_nms_sosh, by = c("x", "y", "id"))
  r_df_nms <- merge(r_df_nms, r_df_nms_comu, by = c("x", "y", "id"))
  r_df_nms <- merge(r_df_nms, r_df_nms_caau, by = c("x", "y", "id"))
  r_df_nms <- merge(r_df_nms, r_df_nms_rhau, by = c("x", "y", "id"))
  
  # Get the bounding box for the current WEA to set individual limits for each panel
  nms_bbox <- st_bbox(NMS[NMS$id == nms, ])
  
  # Store raster data and bounding box together in the list
  raster_list[[nms]] <- list(data = r_df_nms, bbox = nms_bbox)
}

# Combine all raster data into one data frame
r_df_all <- do.call(rbind, lapply(raster_list, function(x) x$data))

# Reshape the data for ggplot (long format)
r_df_long <- tidyr::pivot_longer(r_df_all, cols = c("bfal", "sosh", "comu", "caau", "rhau"), 
                                 names_to = "species", values_to = "value")

r_df_long$species <- factor(r_df_long$species, levels = c("bfal", "sosh", "comu", "caau", "rhau"))
r_df_long$species <- toupper(as.character(r_df_long$species))

# Reorder the factor levels in uppercase
r_df_long$species <- factor(r_df_long$species, levels = c("BFAL", "SOSH", "COMU", "CAAU", "RHAU"))

bfal_q75 <- quantile(bfal_masked, probs = 0.75, na.rm = TRUE)
comu_q75 <- quantile(comu_masked, probs = 0.75, na.rm = TRUE)
caau_q75 <- quantile(caau_masked, probs = 0.75, na.rm = TRUE)
sosh_q75 <- quantile(sosh_masked, probs = 0.75, na.rm = TRUE)
rhau_q75 <- quantile(rhau_masked, probs = 0.75, na.rm = TRUE)

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
# Loop through each WEA and mask the rasters for all species
for (nms in unique(NMS$id)) {
  # Mask the rasters for each species by the current WEA
  bfal_masked_nms <- mask(bfal_masked_fut, NMS_buffered[NMS_buffered$id == nms, ])
  sosh_masked_nms <- mask(sosh_masked_fut, NMS_buffered[NMS_buffered$id == nms, ])
  comu_masked_nms <- mask(comu_masked_fut, NMS_buffered[NMS_buffered$id == nms, ])
  caau_masked_nms <- mask(caau_masked_fut, NMS_buffered[NMS_buffered$id == nms, ])
  rhau_masked_nms <- mask(rhau_masked_fut, NMS_buffered[NMS_buffered$id == nms, ])
  
  # Convert each masked raster to a data.frame
  r_df_nms_bfal <- as.data.frame(bfal_masked_nms, xy = TRUE, na.rm = TRUE)
  r_df_nms_sosh <- as.data.frame(sosh_masked_nms, xy = TRUE, na.rm = TRUE)
  r_df_nms_comu <- as.data.frame(comu_masked_nms, xy = TRUE, na.rm = TRUE)
  r_df_nms_caau <- as.data.frame(caau_masked_nms, xy = TRUE, na.rm = TRUE)
  r_df_nms_rhau <- as.data.frame(rhau_masked_nms, xy = TRUE, na.rm = TRUE)
  
  # Rename columns for each species
  colnames(r_df_nms_bfal)[3] <- "bfal"
  colnames(r_df_nms_sosh)[3] <- "sosh"
  colnames(r_df_nms_comu)[3] <- "comu"
  colnames(r_df_nms_caau)[3] <- "caau"
  colnames(r_df_nms_rhau)[3] <- "rhau"
  
  # Add the 'id' column to each data.frame
  r_df_nms_bfal$id <- nms
  r_df_nms_sosh$id <- nms
  r_df_nms_comu$id <- nms
  r_df_nms_caau$id <- nms
  r_df_nms_rhau$id <- nms
  
  # Combine the data for all  species into one data frame
  r_df_nms <- merge(r_df_nms_bfal, r_df_nms_sosh, by = c("x", "y", "id"))
  r_df_nms <- merge(r_df_nms, r_df_nms_comu, by = c("x", "y", "id"))
  r_df_nms <- merge(r_df_nms, r_df_nms_caau, by = c("x", "y", "id"))
  r_df_nms <- merge(r_df_nms, r_df_nms_rhau, by = c("x", "y", "id"))
  
  # Get the bounding box for the current WEA to set individual limits for each panel
  nms_bbox <- st_bbox(NMS[NMS$id == nms, ])
  
  # Store raster data and bounding box together in the list
  raster_list[[nms]] <- list(data = r_df_nms, bbox = nms_bbox)
}

# Combine all raster data into one data frame
r_df_all_fut <- do.call(rbind, lapply(raster_list, function(x) x$data))

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
r_df_long$value_fut <- r_df_long_fut$value_fut
r_df_long$is_core_fut <- r_df_long_fut$is_core

#subset for plotting
r_df_ch <- subset(r_df_long, id == 'Chumash Heritage')
r_df_mb <- subset(r_df_long, id == 'Monterey Bay')
r_df_cb <- subset(r_df_long, id == 'Cordell Bank')
r_df_gf <- subset(r_df_long, id == 'Greater Farallones')

r_df_ch_core <- subset(r_df_ch, is_core == 1)
r_df_mb_core <- subset(r_df_mb, is_core == 1)
r_df_cb_core <- subset(r_df_cb, is_core == 1)
r_df_gf_core <- subset(r_df_gf, is_core == 1)

r_df_ch_core_fut <- subset(r_df_ch, is_core_fut == 1)
r_df_mb_core_fut <- subset(r_df_mb, is_core_fut == 1)
r_df_cb_core_fut <- subset(r_df_cb, is_core_fut == 1)
r_df_gf_core_fut <- subset(r_df_gf, is_core_fut == 1)

r_df_bfal <- subset(r_df_long, species == 'BFAL')
r_df_sosh <- subset(r_df_long, species == 'SOSH')
r_df_comu <- subset(r_df_long, species == 'COMU')
r_df_caau <- subset(r_df_long, species == 'CAAU')
r_df_rhau <- subset(r_df_long, species == 'RHAU')

mean_labels <- r_df_long %>%
  group_by(id, species) %>%
  summarize(mean_value = mean(value, na.rm = TRUE), .groups = "drop")

# PlOTS ####
a_nms <- ggplot() +
  geom_tile(data = r_df_ch, aes(x = x, y = y, fill = value)) +
  geom_point(data = r_df_ch_core, aes(x = x, y = y), shape = 16, color = "black", size = 1.25, alpha = 1) +
  geom_point(data = r_df_ch_core_fut, aes(x = x, y = y), shape = 16, color = "white", size =1.25, alpha = 1) +
  geom_sf(data = chnms, fill = 'transparent', col = 'black') + 
  facet_wrap(~ species, ncol = 5) +
  scale_fill_viridis_c(name = "HSI", option = "plasma", limits = c(0, 1)) +
  theme_minimal() +
  theme(
    legend.position = 'none',
    strip.text = element_text(size = 10),  # Adjust facet label size if necessary
    axis.text = element_blank()     # Adjust axis label size
  ) +
  labs(
    x = NULL, 
    y = NULL
  )

b_nms <- ggplot() +
  geom_tile(data = r_df_mb, aes(x = x, y = y, fill = value)) +
  geom_point(data = r_df_mb_core, aes(x = x, y = y), shape = 16, color = "black", size = 1.25, alpha = 1) +
  geom_point(data = r_df_mb_core_fut, aes(x = x, y = y), shape = 16, color = "white", size =1.25, alpha = 1) +
  geom_sf(data = mbnms, fill = 'transparent', col = 'black') + 
  facet_wrap(~ species, ncol = 5) +
  guides(fill = 'none') +
  scale_fill_viridis_c(name = "HSI", option = "plasma", limits = c(0, 1)) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10),  # Adjust facet label size if necessary
    axis.text = element_blank()     # Adjust axis label size
  ) +
  labs(
    x = NULL,
    y = NULL
  )

c_nms <- ggplot() +
  geom_tile(data = r_df_cb, aes(x = x, y = y, fill = value)) +
  geom_point(data = r_df_cb_core, aes(x = x, y = y), shape = 16, color = "black", size = 1.25, alpha = 1) +
  geom_point(data = r_df_cb_core_fut, aes(x = x, y = y), shape = 16, color = "white", size =1.25, alpha = 1) +
  geom_sf(data = cbnms, fill = 'transparent', col = 'black') + 
  facet_wrap(~ species, ncol = 5) +
  scale_fill_viridis_c(name = "HSI", option = "plasma", limits = c(0, 1)) +
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

d_nms <- ggplot() +
  geom_tile(data = r_df_gf, aes(x = x, y = y, fill = value)) +
  geom_point(data = r_df_gf_core, aes(x = x, y = y), shape = 16, color = "black", size = 1.25, alpha = 1) +
  geom_point(data = r_df_gf_core_fut, aes(x = x, y = y), shape = 16, color = "white", size =1.25, alpha = 1) +
  geom_sf(data = gfnms, fill = 'transparent', col = 'black') + 
  facet_wrap(~ species, ncol = 5) +
  scale_fill_viridis_c(name = "HSI", option = "plasma", limits = c(0, 1)) +
  theme_minimal() +
  guides() +
  theme_minimal() +
  theme(
    legend.position = 'top',
    strip.text = element_text(size = 10),  # Adjust facet label size if necessary
    axis.text = element_blank()     # Adjust axis label size
  ) +
  labs(
    x = NULL, 
    y = NULL
  )

grid.arrange(d_nms, c_nms, b_nms, a_nms,
             ncol = 1)


# Plot highly suitable habitat richness (# of species with core habitat in a given pixel) 
ch_rich <- r_df_ch %>%
  dplyr::group_by(x, y, id) %>%
  dplyr::summarise(
    richness_hist = sum(as.numeric(as.character(is_core))),       # historical
    richness_fut = sum(as.numeric(as.character(is_core_fut)))     # future
  ) %>%
  dplyr::ungroup()

ch_rich_long <- ch_rich %>%
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

gf_rich <- r_df_gf %>%
  dplyr::group_by(x, y, id) %>%
  dplyr::summarise(
    richness_hist = sum(as.numeric(as.character(is_core))),       # historical
    richness_fut = sum(as.numeric(as.character(is_core_fut)))     # future
  ) %>%
  dplyr::ungroup()

gf_rich_long <- gf_rich %>%
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

richness_nms <- rbind(gf_rich_long, cb_rich_long, mb_rich_long, ch_rich_long)
richness_nms$period <- factor(richness_nms$period, levels = c("Historical", "Future"))

pals <- wesanderson::wes_palette('Zissou1', type = 'continuous')

ch_rich_plot <- ggplot() +
  geom_tile(data = subset(richness_nms, id == 'Chumash Heritage'), aes(x = x, y = y, fill = richness)) +
  geom_sf(data = chnms, fill = 'transparent', lwd = 0.75, col = 'black') +
  facet_wrap(id ~ period, ncol = 2) +
  scale_fill_gradientn(colors = pals, limits = c(0,5)) +
  theme_minimal() +
  guides() +
  theme(
    legend.position = 'none',
    aspect.ratio = 1,
    strip.text = element_blank(),
    plot.margin = margin(0, 0, 0, 0),
    axis.text = element_blank()     # Adjust axis label size
  ) +
  labs(
    x = NULL, 
    y = NULL)

mb_rich_plot <- ggplot() +
  geom_tile(data = subset(richness_nms, id == 'Monterey Bay'), aes(x = x, y = y, fill = richness)) +
  # geom_path(data = mbnms_df, aes(x = X, y = Y, group = L1), color = 'black') +
  geom_sf(data = mbnms, fill = 'transparent', lwd = 0.75, col = 'black') +
  facet_wrap(id ~ period, ncol = 2) +
  scale_fill_gradientn(colors = pals, limits = c(0,5)) +
  theme_minimal() +
  guides() +
  theme(
    legend.position = 'none',
    plot.margin = margin(0, 0, 0, 0),
    aspect.ratio = 1,
    strip.text = element_blank(),
    axis.text = element_blank()     # Adjust axis label size
  ) +
  labs(
    x = NULL, 
    y = NULL)

cb_rich_plot <- ggplot() +
  geom_tile(data = subset(richness_nms, id == 'Cordell Bank'), aes(x = x, y = y, fill = richness)) +
  # geom_path(data = mbnms_df, aes(x = X, y = Y, group = L1), color = 'black') +
  geom_sf(data = cbnms, fill = 'transparent', lwd = 0.75, col = 'black') +
  facet_wrap(id ~ period, ncol = 2) +
  scale_fill_gradientn(colors = pals, limits = c(0,5)) +
  theme_minimal() +
  guides() +
  theme(
    legend.position = 'none',
    aspect.ratio = 1,
    strip.text = element_blank(),
    plot.margin = margin(0, 0, 0, 0),
    axis.text = element_blank()     # Adjust axis label size
  ) +
  labs(
    x = NULL, 
    y = NULL)

gf_rich_plot <- ggplot() +
  geom_tile(data = subset(richness_nms, id == 'Greater Farallones'), aes(x = x, y = y, fill = richness)) +
  # geom_path(data = mbnms_df, aes(x = X, y = Y, group = L1), color = 'black') +
  geom_sf(data = gfnms, fill = 'transparent', lwd = 0.75, col = 'black') +
  facet_wrap(id ~ period, ncol = 2) +
  scale_fill_gradientn(colors = pals, limits = c(0,5)) +
  theme_minimal() +
  guides() +
  theme(
    legend.position = 'none',
    plot.margin = margin(0, 0, 0, 0),
    aspect.ratio = 1,
    strip.text = element_blank(),
    axis.text = element_blank()     # Adjust axis label size
  ) +
  labs(
    x = NULL, 
    y = NULL,
    fill = '# of Species'
  )



# grid.arrange(gf_rich_plot, cb_rich_plot, mb_rich_plot, ch_rich_plot, ncol = 1)


