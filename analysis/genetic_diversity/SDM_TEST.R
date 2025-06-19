library(tidyverse)
library(here)
library(terra)
library(sf)
library(mgcv)
source(here("analysis", "genetic_diversity", "functions_genetic_diversity.R"))
source(here("general_functions.R"))
plotpath <- here("analysis", "genetic_diversity", "plots")

het <- get_het()
coords <- get_coords(sf = TRUE)
hindcast <- rast(here("analysis", "genetic_diversity", "outputs", "hindcast_sdm_sceloporus_occidentalis.tif"))
hindcast_var <- app(hindcast, var, na.rm = TRUE)
# NEED TO CHECK NEW PALEOCLIM DATA FROM CHELSA for cur/lgm
hindcast_mean <- mean(hindcast, na.rm = TRUE) 
stability <- c(hindcast_var, hindcast_mean)
names(stability) <- c("stability_var", "stability_mean")
stability_vals <- extract(stability, coords, ID = FALSE)

mod_df <- read_csv(here("analysis", "genetic_diversity", "outputs", "model_df.csv"))
df <-
  left_join(coords, mod_df) %>%
  bind_cols(stability_vals) %>%
  drop_na(Ho, starts_with("stability"))

pdf(here(plotpath, "hindcast.pdf"), width = 3, height = 3)
gglm("stability_mean", "Ho", df)

gglm("stability_var", "Ho", df)

ggpartial("stability_mean", "Ho", c("stability_mean", "Q"), df)

ggpartial("stability_var", "Ho", c("stability_var", "Q"), df)

ggpartial("csi_past", "Ho", c("csi_past", "Q"), df)
dev.off()

lm(scale(Ho) ~ scale(stability) + scale(Q), data = df) %>% AIC()
lm(scale(Ho) ~ scale(csi_past) + scale(Q), data = df) %>% AIC()
gam(scale(Ho) ~ s(scale(stability), k = 2) + scale(Q), data = df) %>% AIC()
gam(scale(Ho) ~ s(scale(csi_past), k = 2) + scale(Q), data = df) %>% AIC()
