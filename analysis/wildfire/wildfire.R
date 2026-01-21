library(sf)
library(terra)
library(here)
library(tidyverse)
source(here("general_functions.R"))
source(here("analysis", "genetic_diversity", "functions_genetic_diversity.R"))

# Load wildfire data
fire <- rast(here("analysis", "wildfire", "Fire_Summary_Rasters_GeoTiffs", "USGS_Wildland_Fire_Frequency_Raster.tif"))
ca <- get_ca() %>% st_transform(crs(fire))
fire_ca <- crop(fire, ca)

het <- get_het()
nonsyn_het <- read_table(here("analysis", "check_nonsyn", "outputs", "all_nonsynonymous.het")) %>% mutate(Ho_nonsyn = 1 - (`O(HOM)` / `N(NM)`)) %>% select(SampleID = IID, Ho_nonsyn)

coords <- get_coords(sf = TRUE)

model_df <- read_csv(here("analysis", "genetic_diversity", "outputs", "model_df.csv"))
names(model_df)

df <- 
  left_join(coords, het) %>% 
  left_join(nonsyn_het) %>%
  mutate(Ho_resid = resid(lm(Ho_nonsyn ~ Ho))) %>%
  left_join(select(model_df, SampleID, tmean_dif, bio1, csi_past, gHM, NDVI, glacier, Q, fire_severity, fire_frequency, groupveg))

df$recent_fire <- terra::extract(fire, df, ID = FALSE)$USGS_Wildland_Fire_Frequency_Raster

pdf(here("analysis", "wildfire", "wildfire_vs_het.pdf"), width = 4, height = 4)
xlab <- "Burn count\n(number of >10-acre fires in the last century)"
ggplot(df, aes(x = recent_fire, y = Ho)) +
  geom_point(size = 1) +
  geom_smooth(method = "lm") +
  labs(x = xlab, y = "Genome-wide heterozygosity (Ho)") +
  ggpubr::stat_cor() +
  theme_classic() 

ggplot(df, aes(x = log1p(recent_fire), y = Ho)) +
  geom_point(size = 1) +
  geom_smooth(method = "lm") +
  labs(x = paste0("log ", xlab), y = "Genome-wide heterozygosity (Ho)") +
  ggpubr::stat_cor() +
  theme_classic() 
dev.off()

pop_df <- get_pops()

numeric_mod_vars <- c("tmean_dif", "bio1", "csi_past", "Q", "gHM", "NDVI", "Ho", "Ho_nonsyn", "Ho_resid", "recent_fire", "fire_severity", "fire_frequency")
cat_mod_vars <- c("glacier", "groupveg")
unique_df <-
  df %>%
  left_join(pop_df) %>%
  st_transform(3310) %>%
  bind_cols(st_coordinates(df)) %>%
  st_drop_geometry() %>%
  group_by(X, Y) %>%
  mutate(cluster = as.numeric(cluster)) %>%
  summarise(
    across(all_of(c(numeric_mod_vars, "Ho")), ~ mean(.x, na.rm = TRUE)),
    across(all_of(cat_mod_vars), ~ mode_val(.x)),
    .groups = "drop"
  ) %>%
  ungroup() %>%
  # scale
  mutate(across(c(tmean_dif, bio1, csi_past, Q, gHM, NDVI, Ho, Ho_nonsyn, Ho_resid, recent_fire, fire_severity, fire_frequency), scale))

library(nlme)


gls(Ho ~ recent_fire + fire_severity + groupveg + tmean_dif + bio1 + csi_past + Q + gHM + NDVI + glacier, data = unique_df,
    correlation = corExp(form = ~ X + Y, nugget = FALSE)) %>% broom::tidy() %>% mutate(p.value = round(p.value, 3))

gls(Ho ~ recent_fire + fire_severity + groupveg + tmean_dif + bio1 + csi_past + gHM + NDVI + glacier + Q, data = unique_df,
    correlation = corExp(form = ~ X + Y, nugget = FALSE)) %>% broom::tidy() %>% mutate(p.value = round(p.value, 3)) %>% arrange(desc(abs(estimate)))

mod_sf_unique_scaled$fire_recent <- extract(fire_ca, mod_sf_unique_scaled, ID = FALSE)$USGS_Wildland_Fire_Frequency_Raster

gls(Ho ~ bio1 + csi_past + tmean_dif + gHM + NDVI + glacier + fire_severity + 
    groupveg + vdep + Q + fire_recent, data = mod_sf_unique_scaled,
    correlation = corExp(form = ~ x + y, nugget = FALSE)) %>% broom::tidy() %>% filter(p.value < 0.1) %>% mutate(p.value = round(p.value, 3)) %>% arrange(desc(abs(estimate)))

gls(Ho ~ bio1 + csi_past + tmean_dif + gHM + NDVI + glacier + fire_severity + 
    groupveg + vdep + fire_recent, data = mod_sf_unique_scaled,
    correlation = corExp(form = ~ x + y, nugget = FALSE)) %>% broom::tidy() %>% filter(p.value < 0.1) %>% mutate(p.value = round(p.value, 3)) %>% arrange(desc(abs(estimate)))
