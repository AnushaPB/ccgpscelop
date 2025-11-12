library(tidyverse)
library(here)
library(sf)
library(terra)
source(here("general_functions.R"))
source(here("analysis", "genetic_diversity", "functions_genetic_diversity.R"))

het <- get_het()

nonsyn_het <- read_table(here("analysis", "check_nonsyn", "outputs", "all_nonsynonymous.het")) %>% mutate(Ho_nonsyn = 1 - (`O(HOM)` / `N(NM)`)) %>% select(SampleID = IID, Ho_nonsyn)

coords <- get_coords(sf = TRUE)

model_df <- read_csv(here("analysis", "genetic_diversity", "outputs", "model_df.csv"))

df <- 
  left_join(het, coords) %>% 
  left_join(nonsyn_het) %>%
  left_join(dplyr::select(model_df, SampleID, tmean_dif, bio1, csi_past, gHM, NDVI, glacier, Q)) %>%
  mutate(Ho_resid = resid(lm(Ho_nonsyn ~ Ho, data = .))) %>%
  st_as_sf()

df <- left_join(coords, model_df)

unique_df <-
  df %>%
  st_transform(3310) %>%
  bind_cols(st_coordinates(df)) %>%
  st_drop_geometry() %>%
  group_by(X, Y) %>%
  #summarize(across(c(Ho, Ho_nonsyn, Ho_resid, tmean_dif, bio1, csi_past, Q, gHM, NDVI, glacier), ~mean(.x, na.rm = TRUE))) %>%
  summarize(across(c(Ho, tmean_dif, bio1, csi_past, Q, gHM, NDVI, glacier), ~mean(.x, na.rm = TRUE))) %>%
  ungroup() %>%
  # scale
  mutate(across(c(tmean_dif, bio1, csi_past, Q, gHM, NDVI), scale)) %>%
  drop_na()

library(nlme)
gls(Ho ~ tmean_dif + bio1 + csi_past + Q + gHM + NDVI + glacier, data = unique_df,
    correlation = corExp(form = ~ X + Y, nugget = FALSE)) %>% summary()
