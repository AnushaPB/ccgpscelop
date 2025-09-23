
library(tidyverse)
library(here)
library(algatr)
library(terra)
library(sf)
source(here("general_functions.R"))
source(here("analysis", "ibdibe", "functions_ibdibe.R"))

ca_proj <- get_ca() %>% st_transform(3310)
plotpath <- here("analysis", "ibdibe", "plots")

# 1. PREPARE DATA ---------------------------------------------------------------------------

# Read in distance based on SNPs in genes
gendist <- read.csv(here("analysis", "ibdibe", "outputs", "genes_dist.csv"), row.names = 1)
colnames(gendist) <- row.names(gendist) 

# Get coordinates
coords <- get_coords()

# Filter coordinates to those in the gendist matrix and transform to sf object
gendist_coords <- 
  coords %>% 
  filter(SampleID %in% row.names(gendist)) %>%
  mutate(SampleID = factor(SampleID, levels = row.names(gendist))) %>%
  arrange(SampleID) %>% 
  st_as_sf(coords = c("x", "y"), crs = 4326) %>% 
  st_transform(3310)

# Check for disagreement in SampleIDs
setdiff(rownames(gendist), gendist_coords$SampleID)  

# Remove coordinates from gendist
gendist <- gendist[rownames(gendist) %in% gendist_coords$SampleID, colnames(gendist) %in% gendist_coords$SampleID]

# Check that the SampleID in gendist_coords matches the row names of gendist
stopifnot(gendist_coords$SampleID == row.names(gendist))

# Read in environmental data and project to 3310
bioclim <- rast(here("data", "env", "california_chelsa_bioclim_1981-2010_V.2.1.tif"))  %>% project("epsg:3310")
bio1 <- bioclim[["CHELSA_bio1_1981-2010_V.2.1"]]
ndvi <- rast(here("data", "env", "california_ndvi_mean_2000_2020.tif")) %>% project("epsg:3310")

# NDVI has a resolution of ~200 m while bio1 has a resolution of ~1 km, so we resample ndvi to match the resolution of bio1 for the purposes of stacking the environmental layers
ndvi_resampled <- resample(ndvi, bio1, method = "bilinear")
envstack <- c(bio1, ndvi_resampled)
names(envstack) <- c("bio1", "ndvi")


# 2. GDM --------------------------------------------------------------------------------

# Extract environmental data
env <- terra::extract(envstack, gendist_coords, ID = FALSE)

# Run GDM
gdm <- gdm_do_everything(gendist = gendist, coords = gendist_coords, env = env, quiet = TRUE)

# Check coefficients
print(gdm$coeff_df)

# Plot maps (takes a little while)
pdf(here(plotpath, "temp_gdm_outputs.pdf"), width = 5, height = 5)
maps <- gdm_map(gdm$model, envstack,  gendist_coords, plot_vars = TRUE)
dev.off()

# 3. PREDICT GENOMIC OFFSET ---------------------------------------------------------------

# Read in future environmental layers and project to 3310
future <- 
  rast(here("analysis", "gea", "outputs", "scelop_adaptive_env_layers","env_fut_2071-2100_GFDL-ESM4_ssp126_ssp585.tif")) %>%
  project("epsg:3310")

# For this example, use 585 scenario
future <- future[[c("CHELSA_bio1_2071-2100_gfdl-esm4_ssp585_V.2.1", "NDVI")]]
names(future) <- c("bio1", "ndvi")

# Crop future layers
future_resampled <- resample(future, envstack[[1]], method = "bilinear")
future_crop <- crop(future_resampled, envstack)

# Predict genomic offset using GDM
gdm_offset_map <- predict(gdm$model, envstack, time = TRUE, predRasts = future_crop)
range_map <- get_range()
gdm_offset_masked <- mask(gdm_offset_map, range_map)
coord_vals <- terra::extract(gdm_offset_map, gendist_coords, ID = FALSE)[,1]
gdm_offset_masked[gdm_offset_masked > max(coord_vals, na.rm = TRUE)] <- NA  # Set values above max to NA
gdm_offset_masked[gdm_offset_masked < min(coord_vals, na.rm = TRUE)] <- NA  # Set values below min to NA

coords <- 
  gendist_coords %>% 
  mutate(offset = coord_vals)

change_raster <- future_crop[["bio1"]] - envstack[["bio1"]]
residuals_raster <- change_raster
valid_idx <- which(!is.na(terra::values(change_raster)) & !is.na(terra::values(gdm_offset_map)))
resids <- rep(NA, ncell(change_raster))
resids[valid_idx] <- residuals(lm(terra::values(gdm_offset_map)[valid_idx] ~ terra::values(change_raster)[valid_idx]))
values(residuals_raster) <- resids
names(residuals_raster) <- "bio1_change_resids"

coords <- 
  coords %>% 
  mutate(
    bio1_change = terra::extract(change_raster, coords, ID = FALSE)[,1],
    bio1_change_resids = terra::extract(residuals_raster, coords, ID = FALSE)[,1]
  )

# Plot offset
ggraster <- gdm_offset_masked
ggraster[terra::values(ggraster) < min(coords$offset, na.rm = TRUE)] <- NA
ggraster[terra::values(ggraster) > max(coords$offset, na.rm = TRUE)] <- NA
plt1 <- 
  wingen::ggplot_gd(gdm_offset_masked, bkg = ca_proj) +
  geom_sf(data = coords, pch = 21, aes(fill = offset)) +
  scale_fill_gradientn(colors = rev(MetBrewer::met.brewer("Hiroshige", type = "continuous"))) +
  labs(fill = "Genomic\noffset")  +
  theme(
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "mm"),
    legend.position = c(0.9, 0.6),  # (x, y) in npc coords (0-1)
    legend.justification = c(1, 0), 
    legend.key.height = unit(0.5, "cm"),
    legend.key.width = unit(0.5, "cm"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )
  
source(here("analysis", "genetic_diversity", "functions_genetic_diversity.R"))

het <- get_het() %>% right_join(coords, by = "SampleID")
plt2 <-
  ggplot(het, aes(x = offset, y = Ho, col = offset)) +
  geom_point(size = 1.5) +
  geom_smooth(method = "lm", color = "black") +
  labs(x = "Genomic offset", y = "Heterozygosity") +
  theme_classic() +
  theme(legend.position = "none") +
  ggpubr::stat_cor() +
  scale_color_gradientn(colors = rev(MetBrewer::met.brewer("Hiroshige", type = "continuous")))


pdf(here(plotpath, "gdm_offset_map.pdf"), width = 5, height = 8)
cowplot::plot_grid(plt1, plt2, nrow = 1, labels = c("A", "B"))
dev.off()

png(here(plotpath, "gdm_offset_map.png"), width = 8*300, height = 4*300, res = 300)
cowplot::plot_grid(plt1, plt2, nrow = 1, labels = c("A", "B"))
dev.off()

ggraster <- mask(residuals_raster, range_map)
ggraster[terra::values(ggraster) < min(coords$bio1_change_resids, na.rm = TRUE)] <- NA
ggraster[terra::values(ggraster) > max(coords$bio1_change_resids, na.rm = TRUE)] <- NA

domain <- max(abs(c(max(coords$bio1_change_resids, na.rm = TRUE), min(coords$bio1_change_resids, na.rm = TRUE))))

plt_offset <- 
  wingen::ggplot_gd(gdm_offset_masked, bkg = ca_proj) +
  geom_sf(data = coords, pch = 21, aes(fill = offset)) +
  scale_fill_gradientn(colors = rev(MetBrewer::met.brewer("Hiroshige", type = "continuous"))) +
  labs(fill = "Genomic\noffset")  +
  theme(
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "mm"),
    legend.position = c(0.9, 0.6),  # (x, y) in npc coords (0-1)
    legend.justification = c(1, 0), 
    legend.key.height = unit(0.5, "cm"),
    legend.key.width = unit(0.5, "cm"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )

plt_change <- 
  wingen::ggplot_gd(mask(change_raster, range_map), bkg = ca_proj) +
  geom_sf(data = coords, pch = 21, aes(fill = bio1_change)) +
  scale_fill_gradientn(colors = rev(MetBrewer::met.brewer("Hiroshige", type = "continuous"))) +
  labs(fill = "Future\ntemperature\nchange")  +
  theme(
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "mm"),
    legend.position = c(0.9, 0.6),  # (x, y) in npc coords (0-1)
    legend.justification = c(1, 0), 
    legend.key.height = unit(0.5, "cm"),
    legend.key.width = unit(0.5, "cm"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )

plt_resid1 <- 
  wingen::ggplot_gd(ggraster, bkg = ca_proj) +
  geom_sf(data = coords, pch = 21, aes(fill = bio1_change_resids)) +
  scale_fill_gradientn(colors = rev(MetBrewer::met.brewer("Hiroshige", type = "continuous")), limits = c(-domain, domain)) +
  labs(fill = "Residuals")  +
  theme(
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "mm"),
    legend.position = c(0.9, 0.6),  # (x, y) in npc coords (0-1)
    legend.justification = c(1, 0), 
    legend.key.height = unit(0.5, "cm"),
    legend.key.width = unit(0.5, "cm"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )

plt_resid2 <-
  ggplot(coords, aes(x = bio1_change, y = offset, col = bio1_change_resids)) +
  geom_point(size = 1.5) +
  geom_smooth(method = "lm", color = "black") +
  labs(x = "Future temperature change", y = "Genomic offset") +
  theme_classic() +
  theme(legend.position = "none") +
  ggpubr::stat_cor() +
  scale_color_gradientn(colors = rev(MetBrewer::met.brewer("Hiroshige", type = "continuous")), limits = c(-domain, domain))

png(here(plotpath, "gdm_offset_resids_map.png"), width = 8*300, height = 8*300, res = 300)
cowplot::plot_grid(plt_offset, plt_change, plt_resid1, plt_resid2, nrow = 2, labels = c("A", "B", "C", "D"))
dev.off()

source(here("analysis", "genetic_diversity", "functions_genetic_diversity.R"))
het <- 
  get_het() %>% 
  filter(SampleID %in% row.names(gendist)) %>%
  left_join(coords)

plt2 <-
  ggplot(het, aes(x = offset, y = Ho)) +
  geom_point(size = 1.5, aes(col = offset)) +
  geom_smooth(method = "lm", color = "black") +
  labs(x = "Genomic offset", y = "Heterozygosity") +
  theme_classic() +
  theme(legend.position = "none") +
  ggpubr::stat_cor() +
  scale_color_gradientn(colors = rev(MetBrewer::met.brewer("Hiroshige", type = "continuous"))) 
 
png(here(plotpath, "gdm_offset_map.png"), width = 8*300, height = 4*300, res = 300)
cowplot::plot_grid(plt2, plt, nrow = 1, labels = c("C", "D"))
dev.off()

# BACKWARDS OFFSET
tmean_his <- rast(here("data", "env", "tmean_his.tif")) %>% project(envstack)
ndvi_resampled <- resample(ndvi, tmean_his, method = "bilinear")
bwd_stack <- c(tmean_his, ndvi_resampled) 
names(bwd_stack) <- c("bio1", "ndvi")

bwd_crop <- crop(bwd_stack, envstack)
gdm_bwoffset_map <- predict(gdm$model, envstack, time = TRUE, predRasts = bwd_crop)
gdm_bwoffset_masked <- mask(gdm_bwoffset_map, range_map)
coords$bwd_offset <- extract(gdm_bwoffset_masked, gendist_coords, ID = FALSE)[,1]
coords <- left_join(coords, het)

plt <- 
  wingen::ggplot_gd(gdm_bwoffset_masked, bkg = ca_proj) +
  geom_sf(data = coords, pch = 1) +
  scale_fill_gradientn(colors = rev(MetBrewer::met.brewer("Hiroshige", type = "continuous"))) +
  labs(fill = "Genomic\noffset")  +
  theme(
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "mm"),
    legend.position = c(0.9, 0.6),  # (x, y) in npc coords (0-1)
    legend.justification = c(1, 0), 
    legend.key.height = unit(0.5, "cm"),
    legend.key.width = unit(0.5, "cm"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )

plt2 <-
  ggplot(coords, aes(x = bwd_offset, y = Ho)) +
  geom_point(size = 1, aes(col = bwd_offset)) +
  geom_smooth(method = "lm", color = "black") +
  labs(x = "Genomic offset", y = "Heterozygosity") +
  theme_classic() +
  theme(legend.position = "none") +
  ggpubr::stat_cor() +
  scale_color_gradientn(colors = rev(MetBrewer::met.brewer("Hiroshige", type = "continuous"))) 
 
plotpath <- here("analysis", "ibdibe", "plots")
png(here(plotpath, "gdm_backwards_offset_map.png"), width = 8*300, height = 4*300, res = 300)
cowplot::plot_grid(plt2, plt, nrow = 1)
dev.off()
