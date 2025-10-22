library(gradientForest)
library(tidyverse)
library(here)
library(sf)
library(terra)
source("general_functions.R")

# X_now, X_fut: data.frames of predictors (scaled), rows = sites
# Y: data.frame of allele freqs (0..1) per locus at the same sites
# Optionally filter SNPs to candidates
# Read in environmental data and project to 3310
bioclim <- rast(here("data", "env", "california_chelsa_bioclim_1981-2010_V.2.1.tif"))  %>% project("epsg:3310")
bio1 <- bioclim[["CHELSA_bio1_1981-2010_V.2.1"]]
ndvi <- rast(here("data", "env", "california_ndvi_mean_2000_2020.tif")) %>% project("epsg:3310")

# NDVI has a resolution of ~200 m while bio1 has a resolution of ~1 km, so we resample ndvi to match the resolution of bio1 for the purposes of stacking the environmental layers
ndvi_resampled <- resample(ndvi, bio1, method = "bilinear")
envstack <- c(bio1, ndvi_resampled)
names(envstack) <- c("bio1", "ndvi")

dos <- read_table(here("analysis", "gea", "outputs", "nonsyn.raw"))
Y <- dos %>% dplyr::select(starts_with("chr")) 
dos_coords <- get_coords(sf = TRUE) %>% filter(SampleID %in% dos$IID) %>% st_transform(3310)
env <- terra::extract(envstack, dos_coords)

set.seed(1)
gf <- gradientForest(
  data = cbind(Y, X_now),           # GF expects predictors last by default; check args
  predictor.vars = colnames(X_now),
  response.vars  = colnames(Y),
  ntree = 500, maxLevel = 8, corr.threshold = 0.5
)

Xnow_t <- predict(gf, X_now)
Xfut_t <- predict(gf, X_fut)

# Euclidean offset:
offset <- sqrt(rowSums((Xfut_t - Xnow_t)^2))
# Mahalanobis (optional):
S <- cov(Xnow_t)
offset_maha <- sqrt(mahalanobis(Xfut_t, center = Xnow_t, cov = S))
