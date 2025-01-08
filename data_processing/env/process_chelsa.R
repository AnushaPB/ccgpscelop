

library(terra)
library(sf)
library(here)
library(tidyverse)
source(here("general_functions.R"))
ca <- get_ca()

#Name: CHELSA V2
#Source: Karger, D.N., Conrad, O., Böhner, J., Kawohl, T., Kreft, H., Soria-Auza, R.W., Zimmermann, N.E., Linder, P., Kessler, M. (2017): Climatologies at high resolution for the Earth land surface areas. Scientific Data. 4 170122. https://doi.org/10.1038/sdata.2017.122
#Link: https://envicloud.wsl.ch/#/?prefix=chelsa%2Fchelsa_V2%2FGLOBAL%2F
#Downloaded using chelsa/get_chelsa.sh script

#Information from website:
#CHELSA (Climatologies at high resolution for the earth’s land surface areas) is a very high resolution (30 arc sec, ~1km) global downscaled climate data set currently hosted by the Swiss Federal Institute for Forest, Snow and Landscape Research WSL. It is built to provide free access to high resolution climate data for research and application, and is constantly updated and refined.

# Get chelsa files
chelsa_files <- list.files(here("data", "env", "chelsa"), full.names = TRUE)
chelsa <- rast(chelsa_files)

# Crop to California
chelsa_crop <- crop(chelsa, ca)
chelsa_mask <- mask(chelsa_crop, ca)

# Export to a single file
writeRaster(chelsa_mask, here("data", "env", "california_chelsa_bioclim_1981-2010_V.2.1.tif"), overwrite = TRUE)

# Create a PCA raster from the bioclimatic layers:
library(RStoolbox)

# Perform raster PCA on just the bioclimatic data
pca <- rasterPCA(scale(chelsa_mask))

writeRaster(pca$map[[1:3]], here("data", "env", "california_chelsa_bioclim_1981-2010_V.2.1_pca.tif"), overwrite = TRUE)
