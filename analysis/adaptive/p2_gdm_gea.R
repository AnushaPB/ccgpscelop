library(here)
library(algatr)
library(tidyverse)
library(terra)
library(raster)
library(gdm)
library(sf)

source(here("general_functions.R"))
source(here("analysis", "adaptive", "adaptive_index.R"))

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)
path = args[1] # here("analysis", "adaptive", "outputs")
prefix = args[2] # "58-Sceloporus_bio1ndvi_gea_ibs"


# Read in files -----------------------------------------------------------

# Get sampling coordinates and env layers
coords_xy <- get_coords()

ids <- read_tsv(paste0(path, "/", prefix, ".mdist.id"),
                col_names = c("tmp", "SampleID")) %>% dplyr::select(SampleID)
# Check that samples are consistent between dist file and coordinates file
all(coords_xy$SampleID == ids$SampleID) # TRUE

# Retrieve BIO1 and NDVI layers
envlayers <- get_envlayers(future = FALSE)

# Extract and standardize environmental variables and make into dataframe
coords <- get_coords(sf = TRUE)
env <- terra::extract(envlayers$env_pres, coords %>% dplyr::select(geometry))
mod_df <- bind_cols(coords, env)

# =========
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
env <- mod_df_imputed %>% st_drop_geometry() %>% dplyr::select(ID, BIO1, NDVI) %>% as.data.frame()


# Process genetic distances -----------------------------------------------

# Retrieve gendists calculated using only RDA outliers
gendist <- algatr::gen_dist(plink_file = paste0(path, "/", prefix, ".mdist"), 
                            plink_id_file =  paste0(path, "/", prefix, ".mdist.id"), 
                            dist_type = "plink")


# Run GDM -----------------------------------------------------------------

gdm_result <- algatr::gdm_run(
  gendist = as.matrix(gendist),
  coords = coords_xy %>% dplyr::select(x, y) %>% as.matrix(),
  env = env %>% dplyr::select(BIO1, NDVI),
  model = "full",
  scale_gendist = TRUE)

saveRDS(gdm_result, paste0(path, "/GDM_GEA_model.RDS"))