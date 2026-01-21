
library(terra)
library(tidyverse)
library(sf)
library(here)
source(here("general_functions.R"))
path <- here("analysis", "wildfire") # CHANGE

ca_proj_buffer <- get_ca() %>% st_transform(5070) %>% st_buffer(dist = 10000) # buffer to ensure full coverage when cropping rasters

# Historical fire regime data
# Fire Return Interval (FRI)
fri_og <- rast(here(path, "LF2016_FRI_200_CONUS/Tif/LC16_FRI_200.tif"))
fri <- crop(fri_og, ca_proj_buffer)

# Percent Replacement (PFS)
pfs_og <- rast(here(path, "LF2016_PFS_200_CONUS/Tif/LC16_PFS_200.tif"))
pfs <- crop(pfs_og, ca_proj_buffer)

# Vegetation Departure (VDep)
vdep_og <- rast(here(path, "LF2016_VDep_200_CONUS/Tif/LC16_VDep_200.tif"))
vdep <- crop(vdep_og, ca_proj_buffer)

# Activate relevant category for each raster
activeCat(fri)<- "FRI_ALLFIR"
activeCat(vdep) <- "LABEL"

# There are multiple categories for PFS, so create separate rasters for each
pfs_replac <- pfs
pfs_mixed <- pfs
pfs_surfac <- pfs
activeCat(pfs_replac) <- "PRC_REPLAC"
activeCat(pfs_mixed) <- "PRC_MIXED"
activeCat(pfs_surfac) <- "PRC_SURFAC"

# Write out processed rasters
fire_list <- list(
  vdep = vdep,
  fri = fri,
  pfs_replac = pfs_replac,
  pfs_mixed = pfs_mixed,
  pfs_surfac = pfs_surfac
)

walk(names(fire_list), function(x) {
  writeRaster(fire_list[[x]], here(path, paste0("california_", x, ".tif")), overwrite = TRUE)
}, .progress = TRUE)


# Make everything non-NA (set -9999 to NA)
fire_list_final <- map(fire_list, function(x) {
  if (names(x) == "LABEL") {
    message("Transforming VDep raster")

    # Convert non-numeric values to NA
    x_num <- as.numeric(x)

    # Replace -9999 with NA
    x_num[x_num == -9999] <- NA

    return(x_num)
  }

  # Create df from raster
  x_df <- as.data.frame(x, xy = TRUE, na.rm = FALSE)
  x_df$new <- ifelse(x_df[[names(x)]] == -9999, NA, x_df[[names(x)]])

  # Convert from df back to raster
  x_rast <- rast(x_df[, c("x", "y", "new")])
  names(x_rast) <- names(x)

}, .progress = TRUE)


walk(names(fire_list_final), function(x) {
  writeRaster(fire_list_final[[x]], here(path, paste0("california_", x, ".tif")), overwrite = TRUE)
}, .progress = TRUE)

# Aggregate to 1 km resolution
# Create template raster based on ca_proj_buffer at 1 km resolution
template_rast <- rast(ext(ca_proj_buffer), resolution = 1000, crs = crs(ca_proj_buffer))
fire_list_1km <- map(fire_list_final, function(x) {
  resampled <- resample(x, template_rast, method = "bilinear")
  return(resampled)
}, .progress = TRUE)

walk(names(fire_list_1km), function(x) {
  writeRaster(fire_list_1km[[x]], here(path, paste0("california_", x, "_aggregated1km.tif")), overwrite = TRUE)
}, .progress = TRUE)
