
library(tidyverse)
library(here)
library(algatr)
library(terra)
library(sf)
source(here("general_functions.R"))
source(here("analysis", "ibdibe", "functions_gendist.R"))

plotpath <- here("analysis", "ibdibe", "plots")

# 1. PREPARE DATA ---------------------------------------------------------------------------

# Read in distance based on SNPs in genes
gendist <- read.csv(here("analysis", "ibdibe", "outputs", "genes_dist.csv"), row.names = 1)
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

# Check that the SampleID in gendist_coords matches the row names of gendist
stopifnot(gendist_coords$SampleID == row.names(gendist))

# Read in environmental data and project to 3310
bioclim <- rast(here("data", "env", "california_chelsa_bioclim_1981-2010_V.2.1.tif"))  %>% project("epsg:3310")
bio1 <- bioclim[["CHELSA_bio1_1981-2010_V.2.1"]]
ndvi <- rast(here("data", "env", "california_ndvi_mean_2000_2020.tif")) %>% project("epsg:3310")

# NDVI has a resolution of ~200 m while bio1 has a resolution of ~1 km, so we resample bio1 to match the resolution of ndvi for the purposes of stacking the environmental layers
bio1_resampled <- resample(bio1, ndvi, method = "bilinear")
envstack <- c(bio1_resampled, ndvi)
names(envstack) <- c("bio1", "ndvi")


# 2. GDM --------------------------------------------------------------------------------

# Extract environmental data
env <- terra::extract(envstack, gendist_coords, ID = FALSE)

# Run GDM
gdm <- gdm_do_everything(gendist = gendist, coords = gendist_coords, env = env, quiet = TRUE)

# Check coefficients
print(gdm$coeff_df)

# Plot maps (takes a little while)
pdf(here(plotpath, "temp_gdm_outputs.pdf"), width = 5, height = 5)
maps <- gdm_map(gdm$model, envstack,  gendist_coords, plot_vars = TRUE)
dev.off()

# 3. PREDICT GENOMIC OFFSET ---------------------------------------------------------------

# Read in future environmental layers and project to 3310
future <- 
  rast(here("analysis", "gea", "outputs", "scelop_adaptive_env_layers","env_fut_2071-2100_GFDL-ESM4_ssp126_ssp585.tif")) %>%
  project("epsg:3310")

# For this example, use 585 scenario
future <- future[[c("CHELSA_bio1_2071-2100_gfdl-esm4_ssp585_V.2.1", "NDVI")]]
names(future) <- c("bio1", "ndvi")

# Crop future layers
future_resampled <- resample(future, envstack[[1]], method = "bilinear")
future_crop <- crop(future_resampled, envstack)

# Predict genomic offset using GDM
gdm_offset_map <- predict(gdm$model, envstack, time = TRUE, predRasts = future_crop)

# Plot offset
pdf(here(plotpath, "gdm_offset_map.pdf"), width = 5, height = 5)
plot(gdm_offset_map, axes = FALSE, box = FALSE, col = viridis::rocket(100))
dev.off()
