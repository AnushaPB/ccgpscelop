library(tidyverse)
library(sf)
library(here)
source("general_functions.R")
model_df <- read_csv("analysis/genetic_diversity/outputs/model_df.csv")

ca <- get_ca() %>% st_transform(3310)
model_df$resid <- residuals(lm(tmean_dif ~ bio1, data = model_df)) 
model_df <- model_df %>% st_as_sf(coords = c("x", "y"), crs = 4326) %>% st_transform(3310)

plt1 <-
  ggplot(model_df) +
  geom_sf(data = ca) +
  geom_sf(aes(col = bio1)) +
  #scale_color_viridis_c(option = "plasma") +
  scale_color_gradientn(colors = c("#f88901", "#ffb11f", "#ffffbb", "cornflowerblue", "#276ff4")) +
  theme_void() +
  labs(col = "Current temp")

plt2 <-
  ggplot(model_df) +
  geom_sf(data = ca) +
  geom_sf(aes(col = tmean_dif)) +
  scale_color_gradientn(colors = c("#f88901", "#ffb11f", "#ffffbb", "cornflowerblue", "#276ff4")) +
  theme_void() +
  labs(col = "Climate change")

plt3 <-
  ggplot(model_df) +
  geom_sf(data = ca) +
  geom_sf(aes(col = resid)) +
  scale_color_viridis_c(option = "turbo") +
  theme_void() +
  labs(col = "Residuals")


pdf(here("analysis", "genetic_diversity", "climate_change_resids.pdf"), width = 12, height = 4)
cowplot::plot_grid(plt1, plt2, plt3, nrow = 1, align = "hv")
dev.off()
