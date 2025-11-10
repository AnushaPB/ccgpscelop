library(sf)
library(terra)
library(here)
library(tidyverse)
source(here("general_functions.R"))
source(here("analysis", "genetic_diversity", "functions_genetic_diversity.R"))
gdb <- "fire24_1.gdb"

# Read + normalize geometry
fires <- 
  st_read(gdb, layer = "firep24_1", quiet = TRUE) %>%
  # Fix geometry issues
  st_zm(drop = TRUE, what = "ZM") %>%
  st_cast("MULTIPOLYGON", warn = FALSE) %>%
  st_make_valid() %>%
  # Filter to last 100 years
  mutate(FIRE_YEAR = as.integer(`YEAR_`)) %>%
  filter(FIRE_YEAR >= 1920, FIRE_YEAR <= 2020) %>%
  # Transform to a projected CRS for area calculations
  st_transform(3310)                        

# Dissolve by year so a pixel burns at most once per year
# This helps address potential duplicate burn data
fires_year <- fires %>%
  group_by(FIRE_YEAR) %>%
  summarise(do_union = TRUE, .groups = "drop")

# 0/1 field for rasterize
fires_year$one <- 1L

# Terra vector objects
v_year <- vect(fires_year)

# Template raster (1 km)
r <- rast(ext(v_year), resolution = 1000, crs = crs(v_year))

# Rasterize each year to a 0/1 layer, then stack
layers <- map(1:nrow(fires_year), \(i) {
  # max of 0/1 ensures "burned at least once this year" per cell
  rasterize(v_year[i], r, field = "one", fun = "max", background = 0)
}, .progress = TRUE)

# Create stack of all years
s <- rast(layers)
names(s) <- paste0("y", fires_year$FIRE_YEAR)

# Sum across years = fire count per pixel
fire_count <- app(s, sum, na.rm = TRUE)

png(here("analysis", "wildfire",  "fire_count_1km.png"), width = 6, height = 5, units = "in", res = 300)
plot(fire_count, main = "Number of fires per 1 km pixel")
dev.off()

het <- get_het()
nonsyn_het <- read_table(here("analysis", "check_nonsyn", "outputs", "all_nonsynonymous.het")) %>% mutate(Ho_nonsyn = 1 - (`O(HOM)` / `N(NM)`)) %>% select(SampleID = IID, Ho_nonsyn)

coords <- get_coords(sf = TRUE)

df <- left_join(coords, het) %>% left_join(nonsyn_het)
df$fire <- terra::extract(fire_count, df, ID = FALSE)$sum
df$fire[is.na(df$fire)] <- 0
df$Ho_resid <- resid(lm(Ho_nonsyn ~ Ho, data = df))
df$fire_bin <- df$fire > 0

model_df <- read_csv(here("analysis", "genetic_diversity", "outputs", "model_df.csv"))
df <- left_join(df, select(model_df, SampleID, tmean_dif, bio1, csi_past, gHM, NDVI, glacier, Q))

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

ggplot(df, aes(x = fire_bin, y = Ho)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, size = 1) +
  labs(x = xlab) +
  theme_classic()

ggplot(df, aes(x = fire_bin, y = Ho_nonsyn)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, size = 1) +
  labs(x = xlab)+
  theme_classic()

ggplot(df, aes(x = fire_bin, y = Ho_resid)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, size = 1) +
  labs(x = xlab) +
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
  mutate(fire_bin = as.numeric(fire_bin), cluster = as.numeric(cluster)) %>%
  summarize(across(c(Ho, Ho_nonsyn, Ho_resid, fire_bin, fire, tmean_dif, bio1, csi_past, Q, gHM, NDVI, glacier, cluster), mean)) %>%
  mutate(fire_bin = as.factor(round(fire_bin)), cluster = as.factor(round(cluster))) %>%
  ungroup() %>%
  # scale
  mutate(across(c(tmean_dif, bio1, csi_past, Q, gHM, NDVI), scale))

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
