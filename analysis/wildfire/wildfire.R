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

df <- 
  left_join(coords, het) %>% 
  left_join(nonsyn_het) %>%
  mutate(Ho_resid = resid(lm(Ho_nonsyn ~ Ho))) %>%
  left_join(select(model_df, SampleID, tmean_dif, bio1, csi_past, gHM, NDVI, glacier, Q))

df$fire <- terra::extract(fire, df, ID = FALSE)$USGS_Wildland_Fire_Frequency_Raster

pdf(here("analysis", "wildfire", "wildfire_vs_het.pdf"), width = 4, height = 4)
xlab <- "Burn count\n(number of >10-acre fires in the last century)"
ggplot(df, aes(x = fire, y = Ho)) +
  geom_point(size = 1) +
  geom_smooth(method = "lm") +
  labs(x = xlab, y = "Genome-wide heterozygosity (Ho)") +
  ggpubr::stat_cor() +
  theme_classic() 

ggplot(df, aes(x = fire, y = Ho_nonsyn)) +
  geom_point(size = 1) +
  geom_smooth(method = "lm") +
  labs(x = xlab) +
  ggpubr::stat_cor() +
  theme_classic() 

ggplot(df, aes(x = fire, y = Ho_resid)) +
  geom_point(size = 1) +
  geom_smooth(method = "lm") +
  labs(x = xlab) +
  ggpubr::stat_cor() +
  theme_classic() 
dev.off()

pop_df <- get_pops()

unique_df <-
  df %>%
  left_join(pop_df) %>%
  st_transform(3310) %>%
  bind_cols(st_coordinates(df)) %>%
  st_drop_geometry() %>%
  group_by(X, Y) %>%
  mutate(cluster = as.numeric(cluster)) %>%
  summarize(across(c(Ho, Ho_nonsyn, Ho_resid, fire, tmean_dif, bio1, csi_past, Q, gHM, NDVI, glacier, cluster), mean)) %>%
  mutate(cluster = as.factor(round(cluster))) %>%
  ungroup() %>%
  # scale
  mutate(across(c(tmean_dif, bio1, csi_past, Q, gHM, NDVI, Ho, Ho_nonsyn, Ho_resid), scale))

library(nlme)


gls(Ho ~ fire + tmean_dif + bio1 + csi_past + Q + gHM + NDVI + glacier, data = unique_df,
    correlation = corExp(form = ~ X + Y, nugget = FALSE)) %>% summary()

gls(Ho_nonsyn ~ fire + tmean_dif + bio1 + csi_past + Q + gHM + NDVI + glacier, data = unique_df,
    correlation = corExp(form = ~ X + Y, nugget = FALSE)) %>% summary()

gls(Ho_resid ~ fire + tmean_dif + bio1 + csi_past + Q + gHM + NDVI + glacier, data = unique_df,
    correlation = corExp(form = ~ X + Y, nugget = FALSE)) %>% summary()

gls(Ho ~ fire + tmean_dif + bio1 + csi_past + Q + gHM + NDVI + glacier + cluster, data = unique_df,
    correlation = corExp(form = ~ X + Y, nugget = FALSE)) %>% summary()
  
gls(Ho_resid ~ fire + tmean_dif + bio1 + csi_past + Q + gHM + NDVI + glacier + cluster, data = unique_df,
    correlation = corExp(form = ~ X + Y, nugget = FALSE)) %>% summary()

gls(Ho_resid ~ tmean_dif * cluster + fire + bio1 + csi_past + Q + gHM + NDVI + glacier, data = unique_df,
    correlation = corExp(form = ~ X + Y, nugget = FALSE)) %>% summary()
