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
coords <- get_coords(sf = TRUE)

ids <- read_tsv(paste0(path, "/", prefix, ".mdist.id"),
                col_names = c("tmp", "SampleID")) %>% dplyr::select(SampleID)
# Check that samples are consistent between dist file and coordinates file
all(coords_xy$SampleID == ids$SampleID) # TRUE

# IF USING OLD DATASET
# coords_xy <- read_tsv(paste0(path, "/GDM_GEA_coords.txt"))
# coords <- st_as_sf(coords_xy, coords = c("x", "y"), crs = 4326)

# Retrieve BIO1 and NDVI layers
envlayers <- get_envlayers(future = FALSE)


# Process env data --------------------------------------------------------

# Extract environmental variables and make into dataframe
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

gdm_result <- algatr::gdm_run(
  gendist = as.matrix(gendist),
  coords = coords_xy %>% dplyr::select(x, y) %>% as.matrix(),
  env = env %>% dplyr::select(BIO1, NDVI),
  model = "full",
  scale_gendist = TRUE)

saveRDS(gdm_result, paste0(path, "/GDM_GEA_model.RDS"))

# Run GDM only on BIO1 (excluding NDVI)

gdm_bio1 <- algatr::gdm_run(
  gendist = as.matrix(gendist),
  coords = coords_xy %>% dplyr::select(x, y) %>% as.matrix(),
  env = env %>% dplyr::select(BIO1),
  model = "full",
  scale_gendist = TRUE)


# Build plots -------------------------------------------------------------

library(cowplot)
theme_set(theme_cowplot())

### Panels A and B

gdm_model <- gdm_result$model

# Predicted dissimilarity plots
obs <- tidyr::as_tibble(gdm_model$observed) %>% dplyr::rename(observed = value)
pred <- tidyr::as_tibble(gdm_model$predicted) %>% dplyr::rename(predicted = value)
ecol <- tidyr::as_tibble(gdm_model$ecological) %>% dplyr::rename(ecological = value)
dat <- cbind(obs, pred, ecol)
datL <- nrow(dat)
# Get data for overlaid lines
overlayX_ecol <- seq(from = min(dat$ecological), to = max(dat$ecological), length = datL)
overlayY_ecol <- 1 - exp(-overlayX_ecol)
overlayY_pred <- overlayX_pred <- seq(from = min(dat$predicted), to = max(dat$predicted), length = datL)

plot_ecol <-
  ggplot(dat, aes(x = ecological, y = observed)) +
  geom_hex() +
  scale_fill_viridis(option = "D") +
  geom_smooth(method = "lm", color = "red") +
  ggplot2::scale_y_continuous(expand = c(0, 0)) +
  ggplot2::xlab("Predicted ecological distance") +
  ggplot2::ylab("Observed compositional dissimilarity") +
  labs(fill = "Sample\n count") +
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))

plot_pred <-
  ggplot(dat, aes(x = predicted, y = observed)) +
  geom_hex() +
  scale_fill_viridis(option = "D") +
  geom_smooth(method = "lm", color = "red") +
  ggplot2::scale_y_continuous(expand = c(0, 0)) +
  ggplot2::xlab("Predicted compositional distance") +
  ggplot2::ylab("Observed compositional dissimilarity") +
  labs(fill = "Sample\n count") +
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))

a_plot <- cowplot::plot_grid(plot_ecol, plot_pred, nrow = 1, scale = 0.9)
b_plot <- algatr::gdm_plot_isplines(gdm_result$model, scales = "free")

a_b_plot <-
  plot_grid(a_plot, b_plot, nrow = 2, labels = c("A", "B"), vjust = .25)

### Panel C
coords_proj <- get_coords(sf = TRUE) %>%
  st_as_sf(coords = c("x", "y"), crs = 4326) %>%
  st_transform(3310)
envlayers <- get_envlayers(future = FALSE)
env_pres <- envlayers$env_pres
env_proj <- env_pres %>% terra::project("epsg:3310")
env_agg <- mask(env_proj, get_range())
ca_proj <- get_ca() %>% st_transform(3310)

maps <- gdm_map(gdm_result$model, env_agg, coords_proj, plot_vars = TRUE)
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

c_plot <- plot_grid(p_gdm_rainbow, p_vars, nrow = 1, rel_widths = c(2, 1))

full_plot <-
  plot_grid(a_b_plot, c_plot, ncols = 2)
