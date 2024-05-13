library(here)
library(tidyverse)
library(sf)
devtools::load_all()

# check depth and missingness, averaged across all SNPS

# PART 1: CHECKING FULL SNP DATASET (POST DEPTH AND MAF FILTER) -------------------------------------------------
# get data
path <- here("data", "processed_data")
depth <- read_table(here(path, "sample_depth_info.idepth"))
miss <- read_table(here(path, "sample_missing_info.imiss"))

coords <- read_table(here("data", "58-Sceloporus.coords.txt"), col_names = FALSE)
colnames(coords) <- c("INDV", "x", "y")
coords <- sf::st_as_sf(coords, coords = c("x", "y"))

df <-
  left_join(depth, miss) %>%
  right_join(coords)

# plotting depth/missingness
ggplot(df) +
  geom_point(aes(x = F_MISS, y = MEAN_DEPTH)) +
  #geom_vline(xintercept = 0.6, lty = "dashed") +
  theme_classic()

# removing all individuals with greater than 20% missing data
df %>% filter(F_MISS > 0.20) %>% pull(INDV)

# plot where the individuals are
ggplot(df) +
  geom_sf(aes(geometry = geometry)) +
  geom_sf(data = filter(df, F_MISS > 0.20), aes(geometry = geometry), col = "red")
