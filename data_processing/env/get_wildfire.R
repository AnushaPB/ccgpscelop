library(here)
library(terra)
library(tidyverse)
library(sf)
source(here("general_functions.R"))

# Combined wildland fire datasets for the United States and certain territories, 1800s-Present (summary rasters) # Link from: https://www.sciencebase.gov/catalog/item/61aa5483d34eb622f699df85

# URL from ScienceBase (temporary signed S3 link)
url <- "https://prod-is-usgs-sb-prod-content.s3.amazonaws.com/61aa5483d34eb622f699df85/Fire_Summary_Rasters_GeoTiffs.zip?AWSAccessKeyId=AKIAI7K4IX6D4QLARINA&Expires=1762898372&Signature=LjGRFScbgELdnI8q7KppjMWyREo%3D"

# Download file
zip_path <- here("data", "env", "Fire_Summary_Rasters_GeoTiffs.zip")
download.file(url, destfile = zip_path)
unzip(zip_path, exdir = here("data", "env"))

# Delete zipped file
file.remove(zip_path)

# Load file
fire <- rast(here("data", "env", "Fire_Summary_Rasters_GeoTiffs", "USGS_Wildland_Fire_Frequency_Raster.tif"))

# Crop to California
ca <- get_ca() %>% st_transform(crs(fire))
fire_ca <- crop(fire, ca)
writeRaster(fire_ca, here("data", "env", "california_fire_frequency.tif"), overwrite = TRUE)


