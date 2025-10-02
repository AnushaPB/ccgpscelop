library(here)
library(algatr)
library(tidyverse)
library(terra)
library(raster)
library(gdm)
library(sf)
library(algatr)
library(here)
library(tidyverse)
library(tidyterra)
library(RColorBrewer)
library(sf)
library(terra)
library(raster)
library(cowplot)
library(RStoolbox)
library(tigris)
library(viridis)

source(here("general_functions.R"))
source(here("analysis", "adaptive", "adaptive_index.R"))
source(here("analysis", "adaptive", "genomic_offset.R"))

path = here("analysis", "adaptive", "outputs")
# path = here("analysis", "adaptive", "outputs_apr22")
prefix = "58-Sceloporus_bio1ndvi_gea_ibs"


# Read in input files -----------------------------------------------------

# Get sampling coordinates and env layers
coords_xy <- get_coords()
coords <- get_coords(sf = TRUE)

ids <- read_tsv(paste0(path, "/", prefix, ".mdist.id"),
                col_names = c("tmp", "SampleID")) %>% dplyr::select(SampleID)
# Check that samples are consistent between dist file and coordinates file
all(coords_xy$SampleID == ids$SampleID) # TRUE

# IF USING OLD DATASET
# coords_xy <- read_tsv(paste0(path, "/GDM_GEA_coords.txt"))
# coords <- st_as_sf(coords_xy, coords = c("x", "y"), crs = 4326)

# Retrieve BIO1 and NDVI layers
envlayers <- get_envlayers(future = TRUE)

# plot(envlayers$env_pres[[1]])
# plot(envlayers$env_fut[[1]])

# Retrieve all layers
# bioclim <- rast(here("data", "env", "california_chelsa_bioclim_1981-2010_V.2.1.tif"))
# # bio1 <- bioclim[["CHELSA_bio1_1981-2010_V.2.1"]]
# ndvi <- terra::rast(here("data", "env", "california_ndvi_mean_2000_2020.tif"))
# ndvi <- terra::resample(ndvi, bioclim, method = "bilinear")
# env_pres <- c(bioclim, ndvi)
#
# # env_fut <- terra::rast()
#
# list.files(here("data", "env", "future", "bio"))
# here("data", "env", "future", "bio", )
# env_fut <- terra::rast(here("data", "env", "future", "env_fut_2071-2100_GFDL-ESM4_ssp126_ssp585.tif"))
# names(env_fut) <- c("CHELSA_bio1_2071-2100_gfdl-esm4_ssp126_V.2.1", "CHELSA_bio1_2071-2100_gfdl-esm4_ssp585_V.2.1", "NDVI")


# Process env data --------------------------------------------------------

# Extract environmental variables and make into dataframe
env <- terra::extract(envlayers$env_pres, coords %>% dplyr::select(geometry))
mod_df <- bind_cols(coords, env)

# IMPUTATION
# Grab nearest samples' env values if there are NAs
# Count NAs in each column
have_nas <- colSums(is.na(mod_df))[colSums(is.na(mod_df)) > 0]
have_nas

# Which samples have NAs?
mod_df %>% filter(is.na(BIO1)) # Scelocci_IW3247
mod_df %>% filter(is.na(NDVI)) # Scelocci_IW3247, Scelocci_IW3281, Sceocc_HBS142509

# Impute to the value of the closest sample
mod_sf <- mod_df %>% st_as_sf() %>% st_transform(3310)
close_impute <- function(var, mod_sf){
  map_dbl(1:nrow(mod_sf), \(i){
    if (!is.na(st_drop_geometry(mod_sf[i, var]))) return(mod_sf[i, var][[1]])
    # Get the closest sample that is not NA
    nearest <- st_nearest_feature(mod_sf[i,], drop_na(mod_sf[-i,]), 1)
    # Get the value of the closest sample
    return(drop_na(mod_sf)[nearest, var][[1]])
  })
}

vars <- names(have_nas)
names(vars) <- vars
mod_sf_imputed <- map(vars, ~close_impute(.x, mod_sf), .progress = TRUE)
mod_df_imputed <- bind_cols(dplyr::select(mod_df, -all_of(vars)), bind_cols(mod_sf_imputed))

# Confirm no NAs:
colSums(is.na(mod_df_imputed))[colSums(is.na(mod_df_imputed)) > 0]

# Extract imputed env values for all samples
env <- mod_df_imputed %>%
  st_drop_geometry() %>%
  dplyr::select(ID, BIO1, NDVI) %>%
  as.data.frame()


# Process genetic distances -----------------------------------------------

# Retrieve gendists calculated using only RDA outliers
gendist <- algatr::gen_dist(plink_file = paste0(path, "/", prefix, ".mdist"),
                            plink_id_file =  paste0(path, "/", prefix, ".mdist.id"),
                            dist_type = "plink")


# Run GDM -----------------------------------------------------------------

### Run using gdm package
gendist <- scale01(gendist)
site <- 1:nrow(gendist)
gdmGen <- cbind(site, gendist)
env <- terra::extract(envlayers$env_pres, coords_xy %>% dplyr::select(x, y))

gdmPred <- data.frame(
  site = site,
  x = coords_xy$x,
  y = coords_xy$y,
  env)

gdmData <-
  gdm::formatsitepair(bioData = gdmGen,
                      bioFormat = 3,
                      XColumn = "x",
                      YColumn = "y",
                      siteColumn = "site",
                      predData = gdmPred
                      )

cc <- stats::complete.cases(gdmData)
gdmData <- gdmData[cc, ]
gdmData <- gdmData %>% dplyr::select(-c(s1.ID, s2.ID))
gdm_result_geo <- gdm::gdm(gdmData, geo = TRUE)
gdm_result_nogeo <- gdm::gdm(gdmData, geo = FALSE)

### Run using algatr package
algatr_gdmData <- algatr::gdm_format(gendist,
                                     coords_xy %>% dplyr::select(x, y),
                                     env_pres,
                                     scale_gendist = TRUE)

algatr_gdm_result <- algatr::gdm_run(
  gendist = as.matrix(gendist),
  coords = coords_xy %>% dplyr::select(x, y) %>% as.matrix(),
  env = env %>% dplyr::select(BIO1, NDVI),
  model = "full",
  scale_gendist = TRUE)
# saveRDS(gdm_result, paste0(path, "/GDM_GEA_model.RDS"))

# Run GDM only on BIO1 (excluding NDVI)
# gdm_result <- algatr::gdm_run(
#   gendist = as.matrix(gendist),
#   coords = coords_xy %>% dplyr::select(x, y) %>% as.matrix(),
#   env = env %>% dplyr::select(BIO1),
#   model = "full",
#   scale_gendist = TRUE)

# Run GDM on all env layers
# gdm_allenv <- algatr::gdm_run(
#   gendist = as.matrix(gendist),
#   coords = coords_xy %>% dplyr::select(x, y) %>% as.matrix(),
#   env = env %>% dplyr::select(BIO1),
#   model = "full",
#   scale_gendist = TRUE)


# Look at GDM results -----------------------------------------------------

custom_gdm_map(gdm_result, env_pres, coords_xy %>% dplyr::select(x, y), plot_vars = TRUE)
maps_geo <- gdm_map(gdm_result_geo, env_pres, coords_xy %>% dplyr::select(x, y), plot_vars = TRUE)
maps_nogeo <- gdm_map(gdm_result_nogeo, env_pres, coords_xy %>% dplyr::select(x, y), plot_vars = TRUE)


algatr_maps <- gdm_map(algatr_gdm_result$model, env_pres, coords_xy %>% dplyr::select(x, y), plot_vars = TRUE)

plot(maps_geo$pcaRastRGB)
plot(maps_nogeo$pcaRastRGB)
plot(algatr_maps$pcaRastRGB)

p_gdm_rainbow <- ggplot() +
  geom_sf(data = ca_proj, fill = "lightgrey", color = "NA") +
  geom_spatraster_rgb(data = maps$pcaRastRGB, r = 1, g = 2, b = 3) +
  theme_map() +
  geom_sf(data = coords_proj, pch = 21)

# Loadings 'legend'
p_vars <-
  gdm_plot_vars(map$pcaSamp, map$pcaRast, map$pcaRastRGB,
                coords = coords_xy %>% dplyr::select(x, y),
                x = "PC1", y = "PC2",
                scl = 1, display_axes = FALSE)

# Check future env --------------------------------------------------------

# Extract layers
env_pres <- envlayers$env_pres
env_fut_26 <- terra::subset(envlayers$env_fut, c("CHELSA_bio1_2071-2100_gfdl-esm4_ssp126_V.2.1", "NDVI"))
env_fut_85 <- terra::subset(envlayers$env_fut, c("CHELSA_bio1_2071-2100_gfdl-esm4_ssp585_V.2.1", "NDVI"))

# Make env names consistent
names(env_fut_26) <- names(env_pres)
names(env_fut_85) <- names(env_pres)

# Check that layers look good for each RCP
terra::compareGeom(env_pres[[1]],
                   env_fut_26[[1]],
                   lyrs = TRUE,
                   crs = TRUE,
                   ext = TRUE,
                   rowcol = TRUE)
terra::compareGeom(env_pres[[1]],
                   env_fut_85[[1]],
                   lyrs = TRUE,
                   crs = TRUE,
                   ext = TRUE,
                   rowcol = TRUE)


# Run GDM-based offset ----------------------------------------------------

gdmoff_26_geo <- predict_gdm(object = gdm_result_geo,
                         data = env_pres,
                         time = TRUE,
                         predRasts = env_fut_26,
                         filename = here("analysis", "adaptive", "outputs", "58-Sceloporus_GDMoffset_RCP26_IBSdist.tif"),
                         overwrite = TRUE)

gdmoff_26_nogeo <- predict_gdm(object = gdm_result_nogeo,
                             data = env_pres,
                             time = TRUE,
                             predRasts = env_fut_26,
                             filename = here("analysis", "adaptive", "outputs", "58-Sceloporus_GDMoffset_RCP26_IBSdist.tif"),
                             overwrite = TRUE,
                             geo = FALSE)


gdmoff_85_nogeo <- predict_gdm(object = gdm_result_nogeo,
                         data = env_pres,
                         time = TRUE,
                         predRasts = env_fut_26,
                         filename = here("analysis", "adaptive", "outputs", "58-Sceloporus_GDMoffset_RCP85_IBSdist.tif"),
                         overwrite = TRUE,
                         geo = FALSE)



# Plot offset -------------------------------------------------------------

ca <- get_ca()
range <- get_range()
# Mask to range
mask_gdmoff_85 <- terra::mask(gdmoff_85 %>% terra::project(crs(range)), range)
mask_gdmoff_26 <- terra::mask(gdmoff_26 %>% terra::project(crs(range)), range)

p_gdm_offset_26 <- ggplot() +
  geom_sf(data = ca, fill = "lightgrey") +
  geom_spatraster(data = mask_gdmoff_26) +
  # scale_fill_viridis_c(option = "F", na.value = "transparent", direction = -1) +
  scale_fill_distiller(palette = "Spectral", na.value = "transparent", direction = -1) +
  geom_sf(data = ca, fill = NA, size = 0.1) +
  xlab("Longitude") +
  ylab("Latitude") +
  ggplot2::guides(fill = guide_legend(title = "GDM-based genomic offset (RCP 2.6)")) +
  # facet_grid(~lyr) +
  theme_map() +
  theme(panel.grid = element_blank(), plot.background = element_blank(),
        panel.background = element_blank(), strip.text = element_text(size = 11))

# p_gdm_offset_85 <-
  ggplot() +
  geom_sf(data = ca, fill = "lightgrey") +
  geom_spatraster(data = gdmoff_85_nogeo) +
  # scale_fill_viridis_c(option = "F", na.value = "transparent", direction = -1) +
  scale_fill_distiller(palette = "Spectral", na.value = "transparent", direction = -1) +
  geom_sf(data = ca, fill = NA, size = 0.1) +
  xlab("Longitude") +
  ylab("Latitude") +
  ggplot2::guides(fill = guide_legend(title = "GDM-based genomic offset (RCP 8.5)")) +
  # facet_grid(~lyr) +
  theme_map() +
  theme(panel.grid = element_blank(), plot.background = element_blank(),
        panel.background = element_blank(), strip.text = element_text(size = 11))

p_gdm_offset_26
p_gdm_offset_85
# plot_grid(p_gdm_offset_26, p_gdm_offset_85, nrow = 1)


# Offset toy example ------------------------------------------------------

# predDiss <- predict_gdm(gdm_result_geo, gdmData)

##time example
rastFile <- system.file("./extdata/swBioclims.grd", package="gdm")
envRast <- terra::rast(rastFile)

##make some fake climate change data
futRasts <- env_pres
# Increase temp by 25%
futRasts[[1]] <- futRasts[[1]]+10
# Reduce NDVI by 25%
futRasts[[2]] <- futRasts[[2]]+10

# Offset on bio1 only
timePred <- predict_gdm(gdm_bio1$model, env_pres[[1]], time = TRUE, predRasts = env_fut_85[[1]])
terra::plot(timePred)

# Plot difference between BIO1 present and future layers
new_fut <- terra::resample(env_fut_85[[1]], env_p)

diff <- env_pres[[1]] - env_fut_85[[1]]
plot(diff)

diff_ndvi <- env_pres[[2]] - env_fut_85[[2]]
plot(diff_ndvi)


# Run GDM twice -----------------------------------------------------------

env_vals_pres <- terra::extract(envlayers$env_pres, coords %>% dplyr::select(geometry))
env_vals_fut <- terra::extract(env_fut_85, coords %>% dplyr::select(geometry))

gdm_pres <- algatr::gdm_run(
  gendist = as.matrix(gendist),
  coords = coords_xy %>% dplyr::select(x, y) %>% as.matrix(),
  env = env_vals_pres %>% dplyr::select(BIO1, NDVI),
  model = "full",
  scale_gendist = TRUE)
maps_pres <- gdm_map(gdm_pres$model, env_pres, coords_xy %>% dplyr::select(x, y), plot_vars = TRUE)

gdm_fut <- algatr::gdm_run(
  gendist = as.matrix(gendist),
  coords = coords_xy %>% dplyr::select(x, y) %>% as.matrix(),
  env = env_vals_fut %>% dplyr::select(BIO1, NDVI),
  model = "full",
  scale_gendist = TRUE)
maps_fut <- gdm_map(gdm_fut$model, env_fut_85, coords_xy %>% dplyr::select(x, y), plot_vars = TRUE)

diff <- maps_fut$pcaRastRGB - maps_pres$pcaRastRGB
plot(diff)

#' Scale a raster stack from 0 to 255
#'
#' @param s RasterStack
#'
#' @noRd
#' @export
stack_to_rgb <- function(s) {
  stack_list <- as.list(s)
  new_stack <- terra::rast(purrr::map(stack_list, raster_to_rgb))
  return(new_stack)
}


#' Scale raster from 0 to 255
#'
#' @param r SpatRast
#'
#' @noRd
#' @export
raster_to_rgb <- function(r) {
  rmax <- terra::minmax(r)["max", ]
  rmin <- terra::minmax(r)["min", ]
  if ((rmax - rmin) == 0) {
    r[] <- 255
  } else {
    r <- (r - rmin) / (rmax - rmin) * 255
  }
  return(r)
}

diff_rgb <- stack_to_rgb(diff)

ggplot() +
  geom_sf(data = ca, fill = "lightgrey", color = "NA") +
  geom_spatraster_rgb(data = diff_rgb, r = 1, g = 2, b = 3) +
  theme_map()


# Run predict GDM function line by line -----------------------------------

# Specify params
object = gdm_result$model
data = env_pres
time = TRUE
predRasts = env_fut_85
filename = here("analysis", "adaptive", "outputs", "TEST.tif")
overwrite = TRUE

# Check that present and future layers are stackable
terra::compareGeom(data, predRasts, lyrs=TRUE, crs=TRUE, ext=TRUE, rowcol=TRUE)

if (time) {
  for(i in 1:terra::nlyr(data)){
    if(names(data)[i]!=names(predRasts)[i]){
      stop("Layer names do not match the variables used to fit the model.")
    }
  }
  if (object$geo) {
    if(terra::nlyr(data)!=length(object$predictors)-1 | terra::nlyr(predRasts)!=length(object$predictors)-1){
      stop("Number of variables supplied for prediction does not equal the number used to fit the model.")
    }
  } else {
    if(terra::nlyr(data)!=length(object$predictors) | terra::nlyr(predRasts)!=length(object$predictors)){
      stop("Number of variables supplied for prediction does not equal the number used to fit the model.")
    }
  }

  # sets the correct names to the data
  names(data) <- paste0("s1.", names(data))
  names(predRasts) <- paste0("s2.", names(predRasts))

  # stack all the raster layers to for prediction
  data <- c(
    data,
    predRasts
  )
}

# makes the prediction based on the data object
gdm_predict <- function(mod, dat, raster = FALSE, ...) {
  nr <- nrow(dat)
  predicted <- rep(0, times = nr)

  # convert to matrix once
  dat <- as.matrix(dat)
  # if predicting with rasters get xy from interpolate and add constants
  if (raster) {
    const <- matrix(0L, nrow = nr, ncol = 2)
    colnames(const) <- c("distance", "weights")
    xy_cols <- dat[, 1:2]
    colnames(xy_cols) <- c("s1.xCoord", "s1.yCoord")
    dat <- cbind(const, xy_cols, dat)
  }

  z <- .C( "GDM_PredictFromTable",
           dat,
           as.integer(mod$geo),
           as.integer(length(mod$predictors)),
           as.integer(nr),
           as.double(mod$knots),
           as.integer(mod$splines),
           as.double(c(mod$intercept, mod$coefficients)),
           preddata = as.double(predicted),
           PACKAGE = "gdm")

  return(
    z$preddata
  )
}

# if a time prediction, maps the predicted values to a raster and returns
# the layer, otherwise returns a dataframe of the predicted values
if (time) {
  # predict using gdm model and terra package using terra::interpolate to get xy too
  output <- terra::interpolate(
    object = data,
    model = object,
    fun = gdm_predict,
    xyNames = c("s2.xCoord", "s2.yCoord"),
    raster = TRUE,
    na.rm = TRUE,
    filename = filename,
    ...
  )

  return(output)

} else{
  # predict using a data.frame
  output <- gdm_predict(
    mod = object,
    dat = data,
    raster = FALSE
  )

  return(output)
}
}

# is it a raster object
.is_raster <- function(x){
  z <- class(x)
  return(
    z %in% c("SpatRaster", "RasterStack", "RasterLayer", "RasterBrick", "stars")
  )
}

# check for r
.check_rast <- function(r, name = "r"){
  if(!methods::is(r, "SpatRaster")){
    tryCatch(
      {
        r <- terra::rast(r)
      },
      error = function(cond) {
        message(sprintf("'%s' is not convertible to a terra SpatRaster object!", name))
        message(sprintf("'%s' must be a SpatRaster, stars, or Raster* object.", name))
      }
    )
  }


# Look at RDA results -----------------------------------------------------



