library(here)
library(vcfR)
library(tidyverse)
library(algatr)
library(sf)

source(here("general_functions.R"))
source(here("analysis", "adaptive", "adaptive_index.R"))

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)
gea_method = args[1] # "pca" or "bio1ndvi"
nPC = args[2]
path = args[3] # here("analysis", "adaptive", "outputs")
prefix = args[4] # "58-Sceloporus_bio1ndvi_gea" or "58-Sceloporus_bio1ndvi_gea_genes_nonsyn"


# Read in files -----------------------------------------------------------

# Get genetic data and process it
vcf <- vcfR::read.vcfR(paste0(path, "/", prefix, ".vcf"))
gen <- algatr::vcf_to_dosage(vcf)
gen <- algatr::simple_impute(x = gen, FUN = median)

# Get environmental layers
if (gea_method == "pca") envlayers <- import_env_files(type = "rasterPCs", future = FALSE)
if (gea_method == "bio1ndvi") envlayers <- import_env_files(type = "ind_layers", future = FALSE)

# Get sampling coordinates
coords <- get_coords()

# Check ordering and matches
genID <- colnames(vcf@gt[, -1])
overlap <- coords$SampleID %in% genID
# if(!all(overlap)){warning("Missing genetic data for: ", paste(coords$SampleID[!overlap]), ", removing coordinate data for these individuals...")}
coordsF <- coords[overlap,]
coordsF <- coordsF[match(genID, coordsF$SampleID),]
if(!all(coordsF$SampleID == genID)){warning("Order of samples in coordinates and genetic data do not match")}
coords <- coordsF

# Extract and standardize environmental variables and make into dataframe
env <- raster::extract(envlayers$env_pres, coords %>% dplyr::select(x, y))
env <- data.frame(env)

# mod_df <- 
#   get_coords() %>% 
#   left_join(data.frame(env)) %>%
#   # REMOVE CHANNEL ISLAND SAMPLE
#   filter(SampleID != "Scelocci_CCGPMC_MW01-3-14") %>%
#   filter(SampleID != "Scelocci_CHI1382_DAW5-46-21") 


# Scelocci_IW3247 removed because NA env values
if (any(is.na(env))) {
  warning("Missing values found in env data, removing rows with NAs")
  na_env <- env
  gen <- gen[complete.cases(na_env), ]
  coords <- coords[complete.cases(na_env), ]
  # Must come last
  env <- env[complete.cases(na_env), ]
}

# =========
# Grab nearest samples' env values if there are NAs
# Count NAs in each column
have_nas <- colSums(is.na(env))[colSums(is.na(env)) > 0]
have_nas

# Impute to the value of the closest sample
mod_sf <- env %>% st_as_sf() %>% st_transform(3310)
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

env <- scale(env, center = TRUE, scale = TRUE)
# When only one env layer provided, env colnames will be named simply 'env' which is not informative
# if (ncol(env) == 1) colnames(env) <- names(env_pcs)
# colnames(env) <- paste("env_", colnames(env), sep = "")


# Run RDA without PCA corr ------------------------------------------------

moddf <- data.frame(env)
f <- as.formula(paste0("gen ~ ", paste(colnames(env), collapse = "+")))
mod <- vegan::rda(f, data = moddf)

# Run RDA with PCA corr ---------------------------------------------------

pcres <- stats::prcomp(gen)
pc <- pcres$x[, 1:nPC]
moddf_pc <- data.frame(env, pc)
f_pc <- as.formula(paste0("gen ~ ", paste(colnames(env), collapse = "+"), "+ Condition(", paste(colnames(pc), collapse = "+"), ")"))
mod_pc <- vegan::rda(f_pc, data = moddf_pc)

# Export files ------------------------------------------------------------

saveRDS(mod, paste0(path, "/", prefix, "_mod.RDS"))
saveRDS(mod_pc, paste0(path, "/", prefix, "_modPCs.RDS"))

export_rda_files(mod = mod, output_path = path, suffix = paste0(prefix, "_mod"))
export_rda_files(mod = mod_pc, output_path = path, suffix = paste0(prefix, "_modPCs"))

# write_tsv(coords, here(output_path, "RDA_genes_coords.txt"), col_names = TRUE)
