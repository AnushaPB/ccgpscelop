# Combined wildland fire datasets for the United States and certain territories, 1800s-Present (summary rasters)
# Link from: https://www.sciencebase.gov/catalog/item/61aa5483d34eb622f699df85
wget -c -O ../data/env/Fire_Summary_Rasters_GeoTiffs.zip \
'https://prod-is-usgs-sb-prod-content.s3.amazonaws.com/61aa5483d34eb622f699df85/Fire_Summary_Rasters_GeoTiffs.zip?AWSAccessKeyId=AKIAI7K4IX6D4QLARINA&Expires=1762898372&Signature=LjGRFScbgELdnI8q7KppjMWyREo%3D'
unzip ../data/env/Fire_Summary_Rasters_GeoTiffs.zip

# Load wildfire data
Rscript -e '
library(terra)
library(here)
library(sf)
source(here("general_functions.R
fire <- rast(here("analysis", "wildfire", "Fire_Summary_Rasters_GeoTiffs", "USGS_Wildland_Fire_Frequency_Raster.tif"))
ca <- get_ca() %>% st_transform(crs(fire))
fire_ca <- crop(fire, ca)
