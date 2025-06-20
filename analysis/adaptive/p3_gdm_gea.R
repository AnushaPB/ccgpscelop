library(here)
library(algatr)
library(tidyverse)
library(terra)
library(raster)
library(gdm)

source(here("general_functions.R"))
source(here("analysis", "adaptive", "adaptive_index.R"))

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)
gea_method = args[1] # "pca" or "bio1ndvi"
path = args[2] # here("analysis", "adaptive", "outputs")
prefix = args[3] # "58-Sceloporus_bio1ndvi_gdmgea" or "58-Sceloporus_bio1ndvi_gdmgea_genes_nonsyn"


# Read in files -----------------------------------------------------------

# Get sampling coordinates and env layers
coords_xy <- get_coords()

ids <- read_tsv(here("analysis", "adaptive", "outputs", "58-Sceloporus_bio1ndvi_gea.dist.id"),
                col_names = c("tmp", "SampleID")) %>% dplyr::select(SampleID)
coords_xy <- coords_xy %>% filter(SampleID %in% ids$SampleID)
write_tsv(coords_xy, here("analysis", "adaptive", "outputs", "GDM_GEA_coords.txt"))

if (gea_method == "pca") envlayers <- import_env_files(type = "rasterPCs", future = FALSE)
if (gea_method == "bio1ndvi") envlayers <- import_env_files(type = "ind_layers", future = FALSE)

# Extract and standardize environmental variables and make into dataframe
coords <- get_coords(sf = TRUE)
coords <- coords %>% filter(SampleID %in% ids$SampleID)
env <- terra::extract(envlayers$env_pres, coords %>% dplyr::select(geometry))
mod_df <- bind_cols(coords, env)

# =========
# Grab nearest samples' env values if there are NAs
# Count NAs in each column
have_nas <- colSums(is.na(mod_df))[colSums(is.na(mod_df)) > 0]
have_nas

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

# Retrieve gendists calculated using only RDA outliers
gendist <- algatr::gen_dist(plink_file = here("analysis", "adaptive", "outputs", "58-Sceloporus_bio1ndvi_gea.dist"), 
                            plink_id_file = here("analysis", "adaptive", "outputs", "58-Sceloporus_bio1ndvi_gea.dist.id"), 
                            dist_type = "plink")

gdm_wgeo <- gdm_run(
  gendist = as.matrix(gendist),
  coords = coords_xy %>% dplyr::select(x, y) %>% as.matrix(),
  env = env %>% dplyr::select(BIO1, NDVI),
  model = "full",
  scale_gendist = TRUE)
saveRDS(gdm_wgeo, here("analysis", "adaptive", "outputs", "GDM_GEA_model.RDS"))

# gdm_wgeo_best <- gdm_run(
#   gendist = as.matrix(gendist),
#   coords = coords_xy %>% dplyr::select(x, y) %>% as.matrix(),
#   env = env,
#   model = "best",
#   scale_gendist = TRUE
# )
# vars <- gdm::gdm.varImp(gdmData,
#     geo = FALSE,
#     splines = NULL,
#     nPerm = 10,
#     predSelect = TRUE
#   )




# Run GDM without geo dists
# Create GDM formatted data objects
formatted_data <- 
  gdm_format(
    gendist = gendist, 
    coords = coords_xy %>% as.matrix(), 
    env = env,
    scale_gendist = TRUE, 
    geodist_type = "Euclidean", 
    distPreds = NULL, 
    dist_lyr = NULL,
    gdmPred = TRUE,
    gdmGen = TRUE
    )
  
gdmData <- formatted_data$gdmData
gdmPred <- formatted_data$gdmPred
gdmGen <- formatted_data$gdmGen

# Vector of sites (for individual-based sampling, this is just assigning 1 site to each individual)
site <- 1:nrow(gendist)

# RUN GDM -------------------------------------------------------------------------------------------------------
# Remove any remaining incomplete cases
cc <- stats::complete.cases(gdmData)
if (!all(cc)) {
  gdmData <- gdmData[cc, ]
  warning(paste(sum(!cc), "NA values found in gdmData, removing;", sum(cc), "values remain"))
  }
gdm_nogeo <- gdm::gdm(gdmData, geo = TRUE)
