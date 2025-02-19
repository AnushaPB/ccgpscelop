library(here)
library(vcfR)
library(tidyverse)
library(algatr)

source(here("general_functions.R"))
source(here("analysis", "adaptive", "adaptive_index.R"))

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)
gea_method = args[1] # "pca" or "bio1ndvi"
nPC = args[2]

# Read in files -----------------------------------------------------------

# Get genetic data and process it
if (gea_method == "pca") {
  vcf <- vcfR::read.vcfR(here("analysis", "adaptive", "outputs", "58-Sceloporus_pca_gea_genes.vcf"))
  envlayers <- import_env_files(type = "rasterPCs", future = FALSE, model = NULL, years = NULL, RCP = NULL)
}
if (gea_method == "bio1ndvi") {
  vcf <- vcfR::read.vcfR(here("analysis", "adaptive", "outputs", "58-Sceloporus_bio1ndvi_gea_genes.vcf"))
  envlayers <- import_env_files(type = "ind_layers", future = FALSE, model = NULL, years = NULL, RCP = NULL)
}

gen <- algatr::vcf_to_dosage(vcf)
gen <- algatr::simple_impute(x = gen, FUN = median)
coords <- get_coords(sf = FALSE)

# Check ordering and matches
genID <- colnames(vcf@gt[, -1])
overlap <- coords$SampleID %in% genID
# if(!all(overlap)){warning("Missing genetic data for: ", paste(coords$SampleID[!overlap]), ", removing coordinate data for these individuals...")}
coordsF <- coords[overlap,]
coordsF <- coordsF[match(genID, coordsF$SampleID),]
# if(!all(coordsF$SampleID == genID)){warning("Order of samples in coordinates and genetic data do not match")}
coords <- coordsF

# Extract and standardize environmental variables and make into dataframe
env <- raster::extract(envlayers$env_pres, coords %>% dplyr::select(x, y))
env <- scale(env, center = TRUE, scale = TRUE)
env <- data.frame(env)
# When only one env layer provided, env colnames will be named simply 'env' which is not informative
if (ncol(env) == 1) colnames(env) <- names(env_pcs)
colnames(env) <- paste("env_", colnames(env), sep = "")

# Scelocci_IW3247 removed because NA env values
if (any(is.na(env))) {
  warning("Missing values found in env data, removing rows with NAs")
  na_env <- env
  gen <- gen[complete.cases(na_env), ]
  coords <- coords[complete.cases(na_env), ]
  # Must come last
  env <- env[complete.cases(na_env), ]
}

# Run RDA without PCA corr ------------------------------------------------

moddf <- data.frame(env)
f <- as.formula(paste0("gen ~ ", paste(colnames(env), collapse = "+")))
mod <- vegan::rda(f, data = moddf)

# Run RDA with PCA corr ---------------------------------------------------

pcres <- stats::prcomp(gen)
# stats::screeplot(pcres, type = "barplot", npcs = length(pcres$sdev), main = "PCA Eigenvalues")
pc <- pcres$x[, 1:nPC]
moddf_pc <- data.frame(env, pc)
f_pc <- as.formula(paste0("gen ~ ", paste(colnames(env), collapse = "+"), "+ Condition(", paste(colnames(pc), collapse = "+"), ")"))
mod_pc <- vegan::rda(f_pc, data = moddf_pc)

# Export files ------------------------------------------------------------

if (gea_method == "pca") output_path = here("analysis", "adaptive", "outputs", "RDA_PCA")
if (gea_method == "bio1ndvi") output_path = here("analysis", "adaptive", "outputs", "RDA_bio1_ndvi")

saveRDS(mod, here(output_path, "RDA_geagenes_mod.RDS"))
saveRDS(mod_pc, here(output_path, "RDA_geagenes_modPCs.RDS"))

export_rda_files(mod = mod, output_path = output_path, suffix = "mod")
export_rda_files(mod = mod_pc, output_path = output_path, suffix = "modPCs")

write_tsv(coords, here(output_path, "RDA_genes_coords.txt"), col_names = TRUE)
