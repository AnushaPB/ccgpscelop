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
hindcast <- rast(here("data", "hindcast_sdm_sceloporus_occidentalis.tif"))
pdf(here(plotpath, "hindcast_sdm_sceloporus_occidentalis.pdf"), width = 10, height = 10)
plot(hindcast, col = viridis::turbo(100))
plot(mean(hindcast, na.rm = TRUE), col = viridis::turbo(100), main = "Hindcast SDM")
dev.off()

hindcast_var <- app(hindcast, sd, na.rm = TRUE)
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

plot_coords <- function(df, col, fill_name){
  ggplot(df) +
    geom_sf(data = ca) +
    geom_sf(aes_string(col = col)) +
    scale_color_viridis_c(option = "mako") +
    labs(col = fill_name) +
    theme_void() 
}

pdf(here(plotpath, "hindcast.pdf"), width = 8, height = 8)
plt1 <- cowplot::plot_grid(gglm("stability_mean", "Ho", df) + xlab("Hindcasted SDM stability"), gglm("csi_past", "Ho", df), nrow = 1)
plt2 <- cowplot::plot_grid(plot_coords(df, "stability_mean", "SDM"), plot_coords(df, "csi_past", "CSI"), nrow = 1)
cowplot::plot_grid(plt1, plt2, nrow = 2)
dev.off()

lm(scale(Ho) ~ scale(stability) + scale(Q), data = df) %>% AIC()
lm(scale(Ho) ~ scale(csi_past) + scale(Q), data = df) %>% AIC()
gam(scale(Ho) ~ s(scale(stability), k = 2) + scale(Q), data = df) %>% AIC()
gam(scale(Ho) ~ s(scale(csi_past), k = 2) + scale(Q), data = df) %>% AIC()
