library(terra)
library(here)
library(tidyverse)
ghm_global <- rast(here("data", "env", "ghm_global", "gHM", "gHM.tif"))

csi_past <- rast(here("data", "env", "csi", "Layers", "past", "csi_past.tif"))

csi_past_t <- csi_past %>% project(ghm_global)
csi_past_t2 <- resample(csi_past_t, ghm_global, method = "bilinear")

# Calculate raster correlation
ghm_csi_cor <- layerCor(c(ghm_global, csi_past_t), use = "pairwise.complete.obs")
print(ghm_csi_cor)
