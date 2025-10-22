
library(tidyverse)
library(here)
library(algatr)
library(terra)
library(sf)
source(here("general_functions.R"))
source(here("analysis", "ibdibe", "functions_ibdibe.R"))
tp <- function(x){
  png(here("TEMP.png"))
  plot(x, col = viridis::turbo(100))
  dev.off()
}
ca_proj <- st_read(here("data", "ca_state", "CA_State.shp")) %>% st_transform(3310)

plotpath <- here("analysis", "ibdibe", "plots")

# 1. PREPARE DATA ---------------------------------------------------------------------------
# Format genetic distance matrix
format_dist_helper("nonsyn", "nonsyn_dist")

# Read in distance based on SNPs in genes
gendist <- read.csv(here("analysis", "ibdibe", "outputs", "nonsyn_dist.csv"), row.names = 1)
colnames(gendist) <- row.names(gendist) 

# Get coordinates
coords <- get_coords()

# Filter coordinates to those in the gendist matrix and transform to sf object
gendist_coords <- 
  coords %>% 
  filter(SampleID %in% row.names(gendist)) %>%
  mutate(SampleID = factor(SampleID, levels = row.names(gendist))) %>%
  arrange(SampleID) %>% 
  st_as_sf(coords = c("x", "y"), crs = 4326) %>% 
  st_transform(3310)

# Check for disagreement in SampleIDs
setdiff(rownames(gendist), gendist_coords$SampleID)  

# Remove coordinates from gendist
gendist <- gendist[rownames(gendist) %in% gendist_coords$SampleID, colnames(gendist) %in% gendist_coords$SampleID]

# Check that the SampleID in gendist_coords matches the row names of gendist
stopifnot(gendist_coords$SampleID == row.names(gendist))

# Environmental data
ndvi_unprojected <- rast(here("data", "env", "california_ndvi_mean_2000_2020.tif")) 

global_bioclim <- rast(list.files(here("data", "env", "chelsa"), full.names = TRUE))
ca_bioclim <- crop(global_bioclim, ndvi_unprojected)
bioclim <- ca_bioclim %>% project("epsg:3310")
ndvi <- ndvi_unprojected %>% project("epsg:3310")

# Rename bioclim
names(bioclim) <- gsub("_1981-2010_V.2.1", "", gsub("CHELSA_", "", names(bioclim)))

# NDVI has a resolution of ~200 m while bio1 has a resolution of ~1 km, so we resample ndvi to match the resolution of bio1 for the purposes of stacking the environmental layers
ndvi_resampled <- resample(ndvi, bioclim, method = "bilinear")
current <- c(bioclim, ndvi_resampled)
names(current) <- c(names(bioclim), "ndvi")

# Read in future environmental layers and project to 3310
future_paths <- list.files(here("data", "env", "chelsa_future"), full.names = TRUE)
future_paths <- future_paths[grepl("bio", future_paths)]  # only bioclim
global_future <- rast(future_paths)
ca_future <- crop(global_future, ndvi_unprojected) 
future_bioclim <- ca_future %>% project("epsg:3310")

# Rename future layers
names(future_bioclim) <- gsub("_2071-2100_gfdl-esm4_ssp585_V.2.1", "", gsub("CHELSA_", "", names(future_bioclim)))

future <- c(future_bioclim, current[["ndvi"]])
names(future) <- c(names(future_bioclim), "ndvi")

# SUBSET ENVSTACK TO MATCH FUTURE STACK
# Some bioclim variables are missing in the future dataset, so we subset envstack to only those variables that are present in the future dataset
subcurrent <- current[[names(future)]]

# 3. FUNCTIONS TO RUN GDM ------------------------------------------------------------
run_gdm_test <- function(model_current, model_future, prefix){
  # Mask to CA
  model_envstack <- mask(model_current, ca_proj)
  model_future <- mask(model_future, ca_proj)
  
  # Extract environmental data
  env <- terra::extract(model_envstack, gendist_coords, ID = FALSE)
  sum(is.na(env))
  
  # Run GDM
  gdm <- gdm_do_everything(gendist = gendist, coords = gendist_coords, env = env, quiet = TRUE)
  
  # Check coefficients
  print(gdm$coeff_df)
  
  # Plot maps (takes a little while)
  pdf(here(plotpath, paste0(prefix, "_gdm_outputs.pdf")), width = 5, height = 5)
  gdm_map(gdm$model, model_envstack,  gendist_coords, plot_vars = TRUE)
  print(gdm_plot_isplines(gdm$model, scales = "free_x", coords = gendist_coords, env = env))
  dev.off()
  
  # Predict genomic offset using GDM
  gdm_offset_map <- predict(gdm$model, model_envstack, time = TRUE, predRasts = model_future)
  range_map <- get_range()

  gdm_offset_masked <- mask(gdm_offset_map, range_map)
  # coord_vals <- terra::extract(gdm_offset_map, gendist_coords, ID = FALSE)[,1]
  # gdm_offset_masked[gdm_offset_masked > max(coord_vals, na.rm = TRUE)] <- NA  # Set values above max to NA
  # gdm_offset_masked[gdm_offset_masked < min(coord_vals, na.rm = TRUE)] <- NA  # Set values below min to NA

  png(here(plotpath, paste0(prefix, "_gdm_offset_map.png")), width = 6, height = 6, units = "in", res = 300)
  plot(gdm_offset_masked)
  dev.off()

  return(list(gdm = gdm, gdm_offset_map = gdm_offset_masked))
}

source(here("analysis", "ibdibe", "gdm_plot.R"))

#################################
#  TEST: BIO1 + BIO12 + NDVI    #
#################################
model_current <- subcurrent[[c("bio1", "bio12", "ndvi")]]
model_future <- future[[c("bio1", "bio12", "ndvi")]]

result_bio1_bio12_ndvi <- run_gdm_test(model_current, model_future, prefix = "bio1_bio12_ndvi")

#########################
#  TEST: BIO1 + BIO12   #
#########################
model_current <- subcurrent[[c("bio1", "bio12")]]
model_future <- future[[c("bio1", "bio12")]]

result_bio1_bio12<- run_gdm_test(model_current, model_future, prefix = "bio1_bio12")

# ---------------------------------------------------------------------------------------------
# PCA FUNCTION
run_pca <- function(cur, fut){
  # Run a RasterPCA on model envstack
  # 1) Compute scaling from CURRENT layers
  mu  <- global(cur, "mean", na.rm = TRUE)$mean
  sig <- global(fut, "sd",   na.rm = TRUE)$sd

  # 2) Get values (matrix: ncell x nlayers) and fit PCA using current's mu/sig
  #    prcomp will do the scaling internally when we pass center/scale.
  # center: a logical value indicating whether the variables should be
  #         shifted to be zero centered. Alternately, a vector of length
  #         equal the number of columns of ‘x’ can be supplied.  The
  #         value is passed to ‘scale’.

  # scale.: a logical value indicating whether the variables should be
  #         scaled to have unit variance before the analysis takes place.
  #         The default is ‘FALSE’ for consistency with S, but in general
  #         scaling is advisable.  Alternatively, a vector of length
  #         equal the number of columns of ‘x’ can be supplied.  The
  #         value is passed to ‘scale’.
  X <- values(cur, mat = TRUE)
  cc <- complete.cases(X)
  pca <- prcomp(X[cc, , drop = FALSE], center = mu, scale. = sig, retx = FALSE)

  # 3) Project the first 3 PCs for current & future with the SAME PCA object
  model_current <- predict(cur, pca, index = 1:3)
  model_future  <- predict(fut, pca, index = 1:3)

  return(list(model_current = model_current, model_future = model_future))
}

######################
#  TEST: PCA on BIO  #
######################
bio_names <- grepl("bio", names(subcurrent))
pca_result <- run_pca(subcurrent[[bio_names]], future[[bio_names]])

result_pca_bio <- run_gdm_test(pca_result$model_current, pca_result$model_future, prefix = "pca_bio")

############################
#  TEST: PCA on BIO + NDVI #
############################
pca_result <- run_pca(subcurrent, future)

result_pca_bio_ndvi <- run_gdm_test(pca_result$model_current, pca_result$model_future, prefix = "pca_bio_ndvi")

######################
#  TEST: (PCA on BIO) + NDVI  #
######################
bio_names <- grepl("bio", names(subcurrent))
pca_result <- run_pca(subcurrent[[bio_names]], future[[bio_names]])
model_current <- c(pca_result$model_current, subcurrent[["ndvi"]])
model_future <- c(pca_result$model_future, future[["ndvi"]])

result_pca_bio_ndvi_separate <- run_gdm_test(model_current, model_future, prefix = "pca_bio_ndvi_separate")
