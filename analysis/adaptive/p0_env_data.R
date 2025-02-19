library(here)
library(tidyverse)
library(raster)
library(RStoolbox)
library(terra)

RCP = 2.6 # options are 8.5, 7, or 2.6 from CHELSA V2
model = "gfdl-esm4" # options are "gfdl-esm4" (highest priority) or "ipsl-cm6a-lr"
years = "2071-2100" # options are "2041-2070" or "2071-2100"
if (RCP == 8.5) ssp <- "ssp585"
if (RCP == 7) ssp <- "ssp370"
if (RCP == 2.6) ssp <- "ssp126"

# RasterPCA on present env layers to get model ----------------------------

env_pres <- raster::stack(here("data", "env", "california_chelsa_bioclim_1981-2010_V.2.1.tif"))
# Rename layers (names must be consistent with future too)
names(env_pres) <- c("BIO1", "BIO10", "BIO11", "BIO12", "BIO13",
                        "BIO14", "BIO15", "BIO16", "BIO17", "BIO18", 
                        "BIO19", "BIO2", "BIO3", "BIO4", "BIO5",
                        "BIO6", "BIO7", "BIO8", "BIO9")
env_pc_mod <- RStoolbox::rasterPCA(scale(env_pres))
env_pcs <- raster::stack(env_pc_mod$map)

if (!file.exists(here("data", "env", "california_chelsa_bioclim_1981-2010_V.2.1_pca.tif"))) raster::writeRaster(env_pcs[[1:3]], here("data", "env", "california_chelsa_bioclim_1981-2010_V.2.1_pca.tif"))
# env_proj <- raster::projectRaster(env_pcs$map, crs = "+proj=longlat") # 4326

# Future env layers -------------------------------------------------------

if (model == "gfdl-esm4") env_future <- raster::stack(list.files(paste0(here("data", "env", "future", "envicloud/chelsa/chelsa_V2/GLOBAL/climatologies"), "/", years, "/GFDL-ESM4/", ssp, "/bio"), pattern = "_bio", full.names = TRUE))
if (model == "ipsl-cm6a-lr") env_future <- raster::stack(list.files(here("data", "env", "future", "bio"), pattern = ssp, full.names = TRUE))

message("There are ", nlayers(env_pres), " present env layers and ", nlayers(env_future), " future env layers")

# Rename vars to keep consistent with env_pres
names(env_future) <- c("BIO1", "BIO10", "BIO11", "BIO12", "BIO13",
                       "BIO14", "BIO15", "BIO16", "BIO17", "BIO18", 
                       "BIO19", "BIO2", "BIO3", "BIO4", "BIO5",
                       "BIO6", "BIO7", "BIO8", "BIO9")

# Crop and resample to match extent and res of env PCs
cropped <- terra::crop(terra::rast(env_future), terra::ext(env_pc_mod$map[[1]]))
masked <- terra::mask(cropped, env_pc_mod$map[[1]])
resamp <- terra::resample(masked, env_pc_mod$map)

# Predict values of future vars using present rasterPCA -------------------

future_pcs <- terra::predict(masked, model = env_pc_mod$model)
names(future_pcs[[1:3]]) <- c("PC1", "PC2", "PC3")
raster::writeRaster(raster::stack(future_pcs[[1:3]]), paste0(here("data", "env", "future"), "/CHELSA_", years, "_", model, "_", ssp, "_V.2.1_pca.tif"), overwrite = TRUE)

future_pcs_resamp <- terra::predict(resamp, model = env_pc_mod$model)
names(future_pcs_resamp[[1:3]]) <- c("PC1", "PC2", "PC3")
raster::writeRaster(raster::stack(future_pcs_resamp[[1:3]]), paste0(here("data", "env", "future"), "/CHELSA_", years, "_", model, "_", ssp, "_V.2.1_pca_resamp.tif"), overwrite = TRUE)
