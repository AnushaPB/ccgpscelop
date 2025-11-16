library(terra)
library(tidyverse)
library(sf)
library(here)
source(here("general_functions.R"))

coords <- get_coords(sf = TRUE)

frg <- rast("LF2016_FRG_200_CONUS_20250814/LF2016_FRG_200_CONUS/Tif/LC16_FRG_200.tif")
bps <- rast("LF2016_BPS_200_CONUS/Tif/LC16_BPS_200.tif")
vdep <- rast("LF2016_VDep_200_CONUS/Tif/LC16_VDep_200.tif")
fri <- rast("LF2016_FRI_200_CONUS/Tif/LC16_FRI_200.tif")
pfs <- rast("LF2016_PFS_200_CONUS/Tif/LC16_PFS_200.tif")

activeCat(frg) <- "FRG_NEW"
activeCat(vdep) <- "LABEL"
activeCat(bps) <- "GROUPVEG"
activeCat(fri) <- "FRI_ALLFIR"
fri_replace <- fri
fri_mixed <- fri
fri_surfac <- fri
activeCat(fri_replace) <- "FRI_REPLAC"
activeCat(fri_mixed) <- "FRI_MIXED"
activeCat(fri_surfac) <- "FRI_SURFAC"

pfs_replac <- pfs
pfs_mixed <- pfs
pfs_surfac <- pfs
activeCat(pfs_replac) <- "PRC_REPLAC"
activeCat(pfs_mixed) <- "PRC_MIXED"
activeCat(pfs_surfac) <- "PRC_SURFAC"

buffer_n_extract <- function(r, points, buffer_dist = 1000) {
  # Buffer points and extract raster values
  points_buffer <- st_buffer(points, dist = buffer_dist) 
  vals <- terra::extract(r, points_buffer, na.rm = TRUE)
  names(vals) <- c("ID", "value")

  # For each ID (point buffer), get the majority value
  majority_vals <- 
    vals %>%
    as_tibble() %>%
    group_by(ID, value) %>%
    summarise(n = n()) %>%
    filter(value != -9999, value != "Water") %>%
    filter(n == max(n)) %>%
    ungroup() %>%
    select(ID, value) %>%
    # For IDs with multiple majority values, take the numeric one over the character one
    group_by(ID) %>%
    filter(ifelse(n() > 1 & any(!is.na(as.numeric(as.character(value)))), !is.na(as.numeric(as.character(value))), TRUE)) %>%
    slice(1)
  
  # Combine back to points so that buffers that have no value that is not -9999 or water can get NA values
  majority_vals_complete <- 
    points %>% 
    mutate(ID = row_number()) %>%
    left_join(majority_vals, by = "ID")
  
  points[[names(r)[1]]] <- majority_vals_complete$value
  
  return(points)
}

# Buffer coordinates and extract majority value within buffer
frg_df <- buffer_n_extract(frg, coords)
bps_df <- buffer_n_extract(bps, coords)
vdep_df <- buffer_n_extract(vdep, coords)
fri_df <- buffer_n_extract(fri, coords)
fri_replace_df <- buffer_n_extract(fri_replace, coords)
fri_mixed_df <- buffer_n_extract(fri_mixed, coords)
fri_surfac_df <- buffer_n_extract(fri_surfac, coords)
pfs_replac_df <- buffer_n_extract(pfs_replac, coords)
pfs_mixed_df <- buffer_n_extract(pfs_mixed, coords)
pfs_surfac_df <- buffer_n_extract(pfs_surfac, coords)

# Make order list of frq values I-A, I-B, I-C, II-A, II-B, II-C, etc. using purrr
library(purrr)
frg_vals <- purrr::map(c("I", "II", "III", "IV", "V"), ~paste0(.x, c("-A", "-B", "-C"))) %>% unlist()
frg_cols <- viridis::turbo(length(frg_vals))
names(frg_cols) <- frg_vals

fire_df <-
  frg_df  %>%
  rename(frg = FRG_NEW) %>%
  mutate(
    frg = na_if(frg, "-9999"),
    frg = na_if(frg, "-1111"),
    frg = na_if(frg, "11"),
    frg = na_if(frg, "12"),
    frg = na_if(frg, "31"),
    frg_group = case_when(
      str_starts(frg, "I-")   ~ 1,
      str_starts(frg, "II-")  ~ 2,
      str_starts(frg, "III-") ~ 3,
      str_starts(frg, "IV-")  ~ 4,
      str_starts(frg, "V-")   ~ 5,
      TRUE ~ NA_real_
    ),
    frg_severity = case_when(
      str_detect(frg, "-A") ~ "low",
      str_detect(frg, "-B") ~ "mixed",
      str_detect(frg, "-C") ~ "high",
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(frg = factor(frg, levels = frg_vals)) %>%
  st_as_sf() %>%
  left_join(st_drop_geometry(bps_df), by = "SampleID") %>%
  left_join(st_drop_geometry(vdep_df), by = "SampleID") %>%
  left_join(st_drop_geometry(fri_df), by = "SampleID") %>%
  left_join(st_drop_geometry(fri_replace_df), by = "SampleID") %>%
  left_join(st_drop_geometry(fri_mixed_df), by = "SampleID") %>%
  left_join(st_drop_geometry(fri_surfac_df), by = "SampleID") %>%
  left_join(st_drop_geometry(pfs_replac_df), by = "SampleID") %>%
  left_join(st_drop_geometry(pfs_mixed_df), by = "SampleID") %>%
  left_join(st_drop_geometry(pfs_surfac_df), by = "SampleID") %>%
  mutate(across(c(starts_with("FRI_"), starts_with("PRC_")),  ~na_if(as.numeric(.x), -9999))) %>%
  mutate(VDep = as.numeric(as.character(LABEL))) %>%
  mutate(VDep_cat = ifelse(is.na(VDep), as.character(LABEL), NA)) %>%
  mutate(VDep_nona = ifelse(is.na(VDep), 101, VDep)) %>%
  mutate(FRI_ALLFIR_log = log(FRI_ALLFIR + 1)) %>%
  mutate(PRC_REPLAC_log = log(PRC_REPLAC + 1)) 

# Combine PRC into single variable using PCA
pfs_pca <- 
  fire_df %>%
  st_drop_geometry() %>%
  select(PRC_REPLAC_log, FRI_ALLFIR_log) %>%
  drop_na() %>%
  prcomp(center = TRUE, scale. = TRUE) %>%
  broom::augment() %>%
  select(PRC_PC = .fittedPC1) 

fire_df_pc <- fire_df %>% drop_na(PRC_REPLAC_log, FRI_ALLFIR_log) %>% mutate(PRC_PC= pfs_pca$PRC_PC)

ca <- get_ca()
plotpath <- here("analysis/wildfire/plots")

frg_plot <-
  ggplot() +
  geom_sf(data = ca) +
  geom_sf(data = fire_df, aes(color = frg), size = 2) +
  theme_void() +
  labs(color = "FRG") +
  scale_color_manual(values = frg_cols) 

vdep_plot <-
  ggplot() +
  geom_sf(data = ca) +
  geom_sf(data = drop_na(fire_df, VDep_cat), aes(shape = VDep_cat), size = 2, color = "darkgray") +
  geom_sf(data = drop_na(fire_df, VDep), aes(color = VDep), size = 2) +
  theme_void() +
  labs(color = "Vegetation\nDeparture", shape = "Other type") +
  scale_color_viridis_c(option = "magma") +
  scale_shape_manual(values = c(17, 15, 3, 4))

bps_plot <-
  ggplot() +
  geom_sf(data = ca) +
  geom_sf(data = fire_df, aes(color = GROUPVEG), size = 2) +
  theme_void() +
  labs(color = "BPS") 

fri_plot <-
  ggplot() +
  geom_sf(data = ca) +
  geom_sf(data = fire_df, aes(color = FRI_ALLFIR), size = 2) +
  theme_void() +
  labs(color = "Fire Return\nInterval (years)") +
  scale_color_viridis_c(option = "turbo")

fri_replac_plot <-
  ggplot() +
  geom_sf(data = ca) +
  geom_sf(data = fire_df, aes(color = FRI_REPLAC), size = 2) +
  theme_void() +
  labs(color = "FRI Replace\nInterval (years)") +
  scale_color_viridis_c(option = "turbo")

fri_mixed_plot <-
  ggplot() +
  geom_sf(data = ca) +
  geom_sf(data = fire_df, aes(color = FRI_MIXED), size = 2) +
  theme_void() +
  labs(color = "FRI Mixed\nInterval (years)") +
  scale_color_viridis_c(option = "turbo")

fri_surfac_plot <-
  ggplot() +
  geom_sf(data = ca) +
  geom_sf(data = fire_df, aes(color = FRI_SURFAC), size = 2) +
  theme_void() +
  labs(color = "FRI Surface\nInterval (years)") +
  scale_color_viridis_c(option = "turbo")

pfs_replac_plot <-
  ggplot() +
  geom_sf(data = ca) +
  geom_sf(data = fire_df, aes(color = PRC_REPLAC), size = 2) +
  theme_void() +
  labs(color = "Percent Fire\nSeverity - Replace") +
  scale_color_viridis_c(option = "turbo")

pfs_mixed_plot <-
  ggplot() +
  geom_sf(data = ca) +
  geom_sf(data = fire_df, aes(color = PRC_MIXED), size = 2) +
  theme_void() +
  labs(color = "Percent Fire\nSeverity - Mixed") +
  scale_color_viridis_c(option = "turbo") 

pfs_surfac_plot <-
  ggplot() +
  geom_sf(data = ca) +
  geom_sf(data = fire_df, aes(color = PRC_SURFAC), size = 2) +
  theme_void() +
  labs(color = "Percent Fire\nSeverity - Surface") +
  scale_color_viridis_c(option = "turbo")

prc_pc_plot<- 
  ggplot() +
  geom_sf(data = ca) +
  geom_sf(data = fire_df_pc, aes(color = PRC_PC), size = 2) +
  theme_void() +
  labs(color = "PCA of\nPercent Fire\nSeverity") +
  scale_color_viridis_c(option = "turbo")

png(here(plotpath, "landfire_map.png"), width = 15, height = 15, units = "in", res = 300)
cowplot::plot_grid(frg_plot, vdep_plot, bps_plot, fri_plot, fri_replac_plot, fri_mixed_plot, fri_surfac_plot, pfs_replac_plot, pfs_mixed_plot, pfs_surfac_plot, prc_pc_plot, ncol = 3)
dev.off()


# Collapse
model_df <-
  read_csv(here("analysis/genetic_diversity/outputs/model_df.csv")) %>% 
  left_join(fire_df) %>%
  # Replace NAs in PRC values with 0
  mutate(PRC_REPLAC = ifelse(is.na(PRC_REPLAC),  0, PRC_REPLAC)) %>%
  mutate(PRC_SURFAC = ifelse(is.na(PRC_SURFAC), 0, PRC_SURFAC)) %>%
  mutate(PRC_MIXED = ifelse(is.na(PRC_MIXED), 0, PRC_MIXED)) 

# Check for NAs in FRI
stopifnot(!any(is.na(model_df$FRI_ALLFIR)))
stopifnot(!any(is.na(model_df$FRI_ALLFIR_log)))

# Plot histograms of FRI and PRC variables
pdf(here(plotpath, "landfire_histograms.pdf"), width = 8, height = 8)
par(mfrow = c(3,2))
hist(model_df$FRI_ALLFIR, main = "Fire Return Interval", xlab = "Years")
hist(model_df$FRI_ALLFIR_log, main = "Log Fire Return Interval", xlab = "Log(Years + 1)")
hist(model_df$PRC_REPLAC, main = "Percent Fire Severity - Replacement", xlab = "Percent")
hist(model_df$PRC_REPLAC_log, main = "Log Percent Fire Severity - Replacement", xlab = "Log(Percent + 1)")
hist(model_df$PRC_REPLAC, main = "Percent Fire Severity - Surface", xlab = "Percent")
hist(model_df$PRC_REPLAC_log, main = "Log Percent Fire Severity - Surface", xlab = "Log(Percent + 1)")
dev.off()

cor_mat <- 
  model_df %>% 
  st_drop_geometry() %>% 
  select(FRI_ALLFIR, PRC_REPLAC, PRC_MIXED, PRC_SURFAC, VDep_nona, csi_past, tmean_dif, bio1, Q, glacier, NDVI, gHM, FRI_ALLFIR_log, PRC_REPLAC_log) %>%
  cor(use = "pairwise.complete.obs")
cor_mat[upper.tri(cor_mat)] <- NA
# Pull out pairs with r > 0.6
high_corr <- which(abs(cor_mat) > 0.6, arr.ind = TRUE)
tibble(
  var1 = rownames(cor_mat)[high_corr[, 1]],
  var2 = colnames(cor_mat)[high_corr[, 2]],
  corr = cor_mat[high_corr]
) %>%
  filter(var1 != var2) %>%
  distinct()

model_df_nona <- 
  drop_na(model_df, Ho, FRI_ALLFIR_log, PRC_MIXED, PRC_SURFAC, csi_past, tmean_dif, bio1, Q, glacier, NDVI, gHM) %>%
  # Scale everything
  mutate(across(c(csi_past, tmean_dif, bio1, Q, glacier, NDVI, gHM, FRI_ALLFIR_log, PRC_MIXED, PRC_SURFAC, VDep_nona), scale))
model <- lm(Ho ~ csi_past + tmean_dif + bio1 + Q + glacier + NDVI + gHM + FRI_ALLFIR_log + PRC_SURFAC + PRC_MIXED, data = model_df_nona)
summary(model)

# Dredge model
library(MuMIn)
options(na.action = "na.fail")
dredge_res <- dredge(model, rank = "AICc")
# PUll out top models
top_models <- get.models(dredge_res, subset = delta < 2)
# Print model formulas
dredge_res %>% head()
model_summaries <- lapply(top_models, summary)
model_summaries[[1]]

library(nlme)
mod_vars <- c("csi_past", "tmean_dif", "bio1", "Q", "glacier", "NDVI", "gHM", "FRI_ALLFIR_log", "PRC_REPLAC", "PRC_SURFAC", "PRC_MIXED", "VDep_nona", "Ho")
model_df_unique <- 
  model_df %>% 
  left_join(get_coords(sf = TRUE) %>% 
  bind_cols(st_coordinates(.)) %>%
  st_drop_geometry()) %>%
  group_by(X, Y) %>%
  summarize(across(all_of(mod_vars), \(x) mean(x, na.rm = TRUE))) %>%
  ungroup() %>%
  # Scale variables including Ho for interpretability
  mutate(across(all_of(mod_vars), scale))

library(nlme)
library(ggpubr)
og_mod <- 
  gls(Ho ~ csi_past + tmean_dif + bio1 + Q + glacier + NDVI + gHM, data = model_df_unique,
    correlation = corExp(form = ~ X + Y, nugget = FALSE)
  ) 
mod <- 
  gls(Ho ~ csi_past + tmean_dif + bio1 + Q + glacier + NDVI + gHM + FRI_ALLFIR_log + PRC_SURFAC + PRC_MIXED + Vdep_nona, data = model_df_unique,
    correlation = corExp(form = ~ X + Y, nugget = FALSE)
  ) 
mod_sig <- 
  gls(Ho ~ csi_past + tmean_dif + bio1 + Q + FRI_ALLFIR_log + PRC_SURFAC, data = model_df_unique,
    correlation = corExp(form = ~ X + Y, nugget = FALSE)
  ) 
summary(mod_sig)
summary(og_mod)
summary(mod)

pdf(here(plotpath, "landfire_ho.pdf"), width = 5, height = 5)
ggplot(model_df, aes(x = log(FRI_ALLFIR + 1), y = Ho)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal() +
  ggpubr::stat_cor() +
  labs(x = "Fire Return Interval (years)", y = "Observed Heterozygosity (Ho)")

ggplot(model_df, aes(x = log(PRC_REPLAC + 1), y = Ho)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal() +
  ggpubr::stat_cor() +
  labs(x = "Percent Fire Severity - Replacement", y = "Observed Heterozygosity (Ho)")

ggplot(model_df, aes(x = PRC_SURFAC, y = Ho)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal() +
  ggpubr::stat_cor() +
  labs(x = "Percent Fire Severity - Surface", y = "Observed Heterozygosity (Ho)")

ggplot(fire_df_pc, aes(x = PRC_PC, y = Ho)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal() +
  ggpubr::stat_cor() +
  labs(x = "PC", y = "Observed Heterozygosity (Ho)")

ggplot(model_df, aes(x = VDep_nona, y = Ho)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal() +
  ggpubr::stat_cor() +
  labs(x = "Vegetation Departure", y = "Observed Heterozygosity (Ho)")

ggplot(model_df, aes(x = frg, y = Ho)) +
  geom_boxplot() +
  theme_minimal() +
  ggpubr::stat_cor() +
  labs(x = "Fire Regime Group", y = "Observed Heterozygosity (Ho)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(model_df, aes(x = factor(frg_group), y = Ho)) +
  geom_boxplot() +
  theme_minimal() +
  ggpubr::stat_cor() +
  labs(x = "Fire Regime Group", y = "Observed Heterozygosity (Ho)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(model_df, aes(x = as.numeric(as.character(frg_group)), y = Ho)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal() +
  ggpubr::stat_cor() +
  labs(x = "Fire Regime Group", y = "Observed Heterozygosity (Ho)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(model_df, aes(x = factor(frg_severity), y = Ho)) +
  geom_boxplot() +
  theme_minimal() +
  ggpubr::stat_cor() +
  labs(x = "Fire Regime Group", y = "Observed Heterozygosity (Ho)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

dev.off()

frg_key <- 
  read_csv(here("analysis/wildfire/LF2016_FRG_200_CONUS_20250814/LF2016_FRG_200_CONUS/CSV_Data/LF2016_FRG_200.csv")) %>%
  distinct(FRG = FRG_NEW, FRG_DESC)

tb <- 
  fire_df %>% 
  select(SampleID, FRG = frg, GROUPVEG, VDEP = LABEL) %>% 
  st_drop_geometry() %>%
  left_join(frg_key) %>%
  group_by(FRG, GROUPVEG, VDEP) %>%
  count() %>%
  arrange(desc(n))

write_csv(tb, here("analysis/wildfire/outputs/landfire_attributes.csv"))


# Existing vegetation type > vegetation departure