library(terra)
library(tidyverse)
library(here)
plotpath <- here("analysis", "gea", "plots")
# Load and resample rasters
csi_future <- rast(here("data", "env", "csi", "Layers", "future", "csi_future_ssp585.tif"))
csi_past   <- rast(here("data", "env", "csi", "Layers", "past", "csi_past.tif"))

csi_future_resample <- resample(csi_future, csi_past, method = "bilinear")

csi_stack <- c(csi_future_resample, csi_past)

csi_stack_crop <- crop(csi_stack, ext(-180, 180, -60, 90))  # Crop to exclude Antarctica
csi_proj <- csi_stack_crop  %>% project("EPSG:6933")
csi_stack_agg <- aggregate(csi_proj, fact = 10)  # aggregate to speed up

test <- focalPairs(csi_stack_agg, w = 15, "pearson", use = "complete.pairwise.obs", na.rm = TRUE)
test_mask <- mask(test, csi_stack_agg[[1]])
pdf(here(plotpath, "csi_global.pdf"), width=7, height=5)
plot(csi_stack_agg)
plot(test_mask)
dev.off()



# Define window sizes
window_sizes <- c(3)
cor_rasters <- purrr::map(window_sizes, ~focalPairs(csi_stack, w = .x, "pearson", use = "complete.obs"), .progress = TRUE)
 
# Stack and plot
walk(cor_rasters, plot)

library(sf)

# Load and combine shapefiles
reptiles <- rbind(
  st_read(here("data", "IUCN", "SCALED_REPTILES", "SCALED_REPTILES_PART1.shp"), quiet = TRUE),
  st_read(here("data", "IUCN", "SCALED_REPTILES", "SCALED_REPTILES_PART2.shp"), quiet = TRUE)
)

mammals <- st_read(here("data", "IUCN", "MAMMALS", "MAMMALS_TERRESTRIAL_ONLY.shp"), quiet = TRUE)

# Crop raster to bounding box of all ranges to save time
r_reptiles <- 1 - crop(csi_stack, st_bbox(reptiles))
r_mammals <- 1 - crop(csi_stack, st_bbox(mammals))

# Function to compute Pearson correlation between the 2 raster layers within a polygon
range_cor <- function(i, ranges, r) {
  x <- ranges[i, ]
  vals <- terra::extract(r, x, ID = FALSE)
  if (nrow(vals) < 10 | all(is.na(rowSums(vals)))) return(NA_real_)
  cor(vals[,1], vals[,2], use = "complete.obs")
}

reptile_cors <- map_dbl(1:nrow(reptiles), ~range_cor(.x, ranges = reptiles, r = r_reptiles), .progress = TRUE)
outpath = here("analysis", "gea", "outputs")
writeLines(as.character(reptile_cors), here(outpath, "reptile_csi_cor.txt"))

reptile_df <- reptiles %>% mutate(r = reptile_cors) %>% mutate(group = "Reptiles") %>% st_drop_geometry() %>% dplyr::select(r, group)

mammal_cors <- map_dbl(1:nrow(mammals), ~range_cor(.x, ranges = mammals, r = r_mammals), .progress = TRUE)
writeLines(as.character(reptile_cors), here(outpath, "mammal_csi_cor.txt"))
mammal_df <- mammals %>% mutate(r = mammal_cors) %>% mutate(group = "Mammals") %>% st_drop_geometry() %>% dplyr::select(r, group)

df <- bind_rows(reptile_df, mammal_df)

# Define gradient band breaks and colors
bands <- tibble(
  xmin = c(-Inf, -0.6, -0.3, 0.3, 0.6),
  xmax = c(-0.6, -0.3, 0.3, 0.6, Inf),
  label = c("strong neg", "neg", "neutral", "pos", "strong pos"),
  fill = c("#008080", "#66b2b2", "#ffffff", "#fdd9b5", "#f47c3c")  # teal to white to orange
)
labels <- tibble(
  x = c(-0.8, -0.45, 0.45, 0.8),
  label = c("strongly\nnegative", "weakly\nnegative", "weakly\npositive", "strongly\npositive")
)

mean(mammal_df$r > 0.3, na.rm = TRUE) 
mean(reptile_df$r > 0.3, na.rm = TRUE) 
mean(df$r > 0.3, na.rm = TRUE)

plt1 <-
  ggplot() +
  geom_vline(aes(xintercept = 0), linetype = "dashed") +
  geom_rect(data = bands, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
            inherit.aes = FALSE, alpha = 0.5) +
  scale_fill_identity() +
  geom_text(data = labels, aes(x = x, y = Inf, label = label),
            vjust = 1.5, size = 3, fontface = "italic", inherit.aes = FALSE) +
  geom_histogram(data = reptile_df, aes(x = r), bins = 20, fill = NA, color = "black") +
  theme_classic() +
  labs(x = "Correlation between past and future climate stability", y = "Species count") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1500)) +
  scale_x_continuous(expand = c(0, 0), limits = c(-1, 1)) +
  ggtitle("A. Reptiles")

plt2 <-
  ggplot() +
  geom_vline(aes(xintercept = 0), linetype = "dashed") +
  geom_rect(data = bands, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
            inherit.aes = FALSE, alpha = 0.5) +
  scale_fill_identity() +
  geom_text(data = labels, aes(x = x, y = Inf, label = label),
            vjust = 1.5, size = 3, fontface = "italic", inherit.aes = FALSE) +
  geom_histogram(data = mammal_df, aes(x = r), bins = 20, fill = NA, color = "black") +
  theme_classic() +
  labs(x = "Correlation between past and future climate stability", y = "Species count") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1500)) +
  scale_x_continuous(expand = c(0, 0), limits = c(-1, 1)) +
  ggtitle("A. Mammals")

pdf(here(plotpath, "csi_cor_ranges.pdf"), width = 5, height = 6)
cowplot::plot_grid(plt1, plt2, nrow = 2)
dev.off()

pdf(here(plotpath, "csi_reptile_cor_ranges.pdf"), width = 5, height = 3.5)
plt1
dev.off()

pdf(here(plotpath, "csi_cor_ranges.pdf"), width = 6, height = 3)
 ggplot() +
  geom_vline(aes(xintercept = 0), linetype = "dashed") +
  geom_rect(data = bands, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
            inherit.aes = FALSE, alpha = 0.5) +
  scale_fill_identity() +
  geom_text(data = labels, aes(x = x, y = Inf, label = label),
            vjust = 1.5, size = 3, fontface = "italic", inherit.aes = FALSE) +
  geom_histogram(data = df, aes(x = r, col = group), bins = 20, fill = NA) +
  theme_classic() +
  scale_color_manual(values = c("Reptiles" = "navy", "Mammals" = "blue")) +
  labs(x = "Correlation between past and future climate stability", y = "Species count", col = "") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 3100)) +
  scale_x_continuous(expand = c(0, 0), limits = c(-1, 1))
dev.off()
