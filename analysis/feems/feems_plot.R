
library("here")
library("sf")
library("terra")
library("tidyverse")
source(here("general_functions.R"))

# Load data
ca <- get_ca()
feems <- get_feems(crop = TRUE, nodes = TRUE)
nodes <- feems$nodes
edges <- feems$edges

# Set-up color palette
d <- max(c(abs(min(edges$log_weight, na.rm = TRUE)), abs(max(edges$log_weight, na.rm = TRUE))))
color_palette <- c("#994000", "#CC5800", "#FF8F33", "#FFAD66", "#FFCA99", "#FFE6CC", "#FBFBFB", "#CCFDFF", "#99F8FF", "#66F0FF", "#33E4FF", "#00AACC", "#007A99")

ggplot() +
  geom_sf(data = ca, fill = "#f0f0f0", color="black")   +
  geom_sf(data = edges, aes(col = log_weight), lwd = 1) +
  geom_sf(data = ca, fill = NA, color="black")   +
  geom_sf(data = nodes, aes(cex = size), fill = NA, col = "black", pch = 21) +
  labs(col = "log10(w)", cex = "Sample size") +
  scale_color_gradientn(colors = color_palette, limits = c(-d, d)) +
  theme_void()


# Create an empty raster
feems_vect <- vect(edges)
# Create an empty raster
r <- rast(ext(feems), resolution = 0.01)
# Rasterize the shapefile
r <- rasterize(feems, r, field = "weight")
# Define the moving window 
window <- matrix(1, nrow = 15, ncol = 15)
# Run the moving window analysis
r_moving_window <- focal(r, w = window, fun = mean, na.rm = TRUE)
r_df <- as.data.frame(r_moving_window, xy = TRUE, ID = FALSE)

ggplot() +
  geom_sf(data = ca, fill = NA, color="black")   +
  geom_raster(data = r_df, aes(x = x, y = y, fill = log(focal_mean)), lwd = 1) +
  #geom_sf(data = coords) +
  scale_fill_gradientn(colors = color_palette, limits = c(-d, d), na.value = NA) +
  theme_void()


