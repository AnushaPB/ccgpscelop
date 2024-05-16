library(here)
library(tidyverse)
library(sf)
devtools::load_all()
source(here("general_functions.R"))

path <- here("data", "processed_data")
coords <- get_coords() %>% rename(INDV = SampleID)

# check depth and missingness, averaged across all SNPS

# PART 0: CHECKING RAW DATASET ----------------------------------------------------------------------------------
# get data
depth <- read_table(here(path, "sample_depth_info_raw.idepth"))
miss <- read_table(here(path, "sample_missing_info_raw.imiss"))

df <-
  left_join(depth, miss) %>%
  right_join(coords)

# plotting depth/missingness
ggplot(df) +
  geom_point(aes(x = F_MISS, y = MEAN_DEPTH)) +
  theme_classic()

# PART 1: CHECKING FULL SNP DATASET (POST DEPTH AND MAF FILTER) -------------------------------------------------
# get data
path <- here("data", "processed_data")
depth <- read_table(here(path, "sample_depth_info.idepth"))
miss <- read_table(here(path, "sample_missing_info.imiss"))

df <-
  left_join(depth, miss) %>%
  right_join(coords)

# plotting depth/missingness
ggplot(df) +
  geom_point(aes(x = F_MISS, y = MEAN_DEPTH)) +
  #geom_vline(xintercept = 0.6, lty = "dashed") +
  theme_classic()

# removing all individuals with greater than 60% missing data
df %>% filter(F_MISS > 0.60) %>% nrow()
df %>% filter(F_MISS <= 0.60) %>% nrow()

# plot where the individuals are
ggplot(df) +
  geom_sf(aes(geometry = geometry)) +
  geom_sf(data = filter(df, F_MISS > 0.20), aes(geometry = geometry), col = "red")
