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
coords <- get_coords()

ids <- read_tsv(here("analysis", "adaptive", "outputs", "58-Sceloporus_bio1ndvi_gea.dist.id"),
                col_names = c("tmp", "SampleID")) %>% dplyr::select(SampleID)
coords <- coords %>% filter(SampleID %in% ids$SampleID)
write_tsv(coords, here("analysis", "adaptive", "outputs", "GDM_GEA_coords.txt"))

if (gea_method == "pca") envlayers <- import_env_files(type = "rasterPCs", future = FALSE)
if (gea_method == "bio1ndvi") envlayers <- import_env_files(type = "ind_layers", future = FALSE)
env <- terra::extract(envlayers$env_pres, coords %>% dplyr::select(x, y))

# Retrieve gendists calculated using only RDA outliers
gendist <- algatr::gen_dist(plink_file = here("analysis", "adaptive", "outputs", "58-Sceloporus_bio1ndvi_gea.dist"), 
                            plink_id_file = here("analysis", "adaptive", "outputs", "58-Sceloporus_bio1ndvi_gea.dist.id"), 
                            dist_type = "plink")

# gdm_full <- gdm_run(
#   gendist = as.matrix(gendist),
#   coords = coords %>% dplyr::select(x, y) %>% as.matrix(),
#   env = env,
#   model = "full",
#   scale_gendist = TRUE
# )


# FORMAT DATA ---------------------------------------------------------------------------------------------------
  
  # Create GDM formatted data objects
  formatted_data <- 
    gdm_format(
      gendist = gendist, 
      coords = coords %>% dplyr::select(x, y), 
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
gdm_full <- gdm::gdm(gdmData, geo = FALSE)  

saveRDS(gdm_full, here("analysis", "adaptive", "outputs", "GDM_GEA_model.RDS"))
