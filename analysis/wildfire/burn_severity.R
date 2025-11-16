library(tidyverse)
library(here)
library(terra)
library(sf)
source(here("general_functions.R"))

plotpath <- here("analysis", "wildfire", "plots")

files <- list.files(here("data", "env", "burn_severity"), pattern = ".tif$", full.names = TRUE)
rasters <- lapply(files, rast)

coords <- get_coords(sf = TRUE)

extract_list <- lapply(rasters, function(r) {
  terra::extract(r, coords)
})

band_names <- names(rasters[[1]])  # e.g. "Severity_1984", ...

n_pts <- nrow(coords)

# initialize output matrix for N points × N bands
res_mat <- matrix(NA_real_, nrow = n_pts, ncol = length(band_names))
colnames(res_mat) <- band_names

# Loop through tile extracts and fill in missing values
for (df in extract_list) {
  
  vals <- as.matrix(df[, -1, drop = FALSE])  # drop ID column
  
  # fill NA spots only
  fill_mask <- is.na(res_mat) & !is.na(vals)
  res_mat[fill_mask] <- vals[fill_mask]
}

burn_df <- as_tibble(res_mat) %>%
  mutate(point_id = row_number())

coords_df <- coords %>%
  mutate(point_id = row_number())

# Severity class table
# 0 - Background
# 1	- Unburned to Low
# 2	- Low
# 3	- Moderate
# 4	- High
# 5	- Increased Greenness
# 6	- Non-Mapping Area
coords_burn <- 
  left_join(coords_df, burn_df, by = "point_id") %>%
  # Recode burn severity
  mutate(across(starts_with("Severity_"), ~{
    case_when(
      .x == 1 ~ 0,   # unburned/regrowth -> 0
      .x == 2 ~ 1,   # low -> 1
      .x == 3 ~ 1,   # moderate -> 1
      .x == 4 ~ 2,   # high -> 2
      TRUE ~ NA
    )
  })) %>%
  # Sum all columns to get total burn severity
  rowwise() %>%
  mutate(
    burn_severity_sum = sum(c_across(starts_with("Severity_")), na.rm = TRUE), 
    # Sum binary yes no
    burn_count = sum(!is.na(c_across(starts_with("Severity_")))),
    # Average burn severity
    burn_severity_avg = mean(c_across(starts_with("Severity_")), na.rm = TRUE)
  ) %>%
  ungroup()

png(here(plotpath, "burn_severity_sum.png"), width = 6, height = 6, units = "in", res = 300)
plot(coords_burn["burn_severity_sum"], main = "Cumulative Burn Severity Score", pch = 16)
plot(coords_burn["burn_count"], main = "Cumulative Burn Severity Score", pch = 16)
plot(coords_burn["burn_severity_avg"], main = "Cumulative Burn Severity Score", pch = 16)
dev.off()

source(here("analysis", "genetic_diversity", "functions_genetic_diversity.R"))
model_df <- read_csv(here("analysis", "genetic_diversity", "outputs", "model_df.csv"))
fire_df <- 
  coords_burn %>%
  left_join(model_df, by = "SampleID")

lm(Ho ~ burn_severity_sum + csi_past + tmean_dif + bio1 + Q, data = fire_df) %>% summary()
lm(Ho ~ burn_severity_avg + csi_past + tmean_dif + bio1 + Q, data = fire_df) %>% summary()
lm(Ho ~ burn_count + csi_past + tmean_dif + bio1 + Q, data = fire_df) %>% summary()

# # Function to recode MTBS burn severity
# | Code  | Class Name                                            | Meaning                                                                                                    |
# | ----- | ----------------------------------------------------- | ---------------------------------------------------------------------------------------------------------- |
# | **0** | **Unburned / Increased Growth**                       | No detectable change or vegetation is greener post-fire (regrowth, seasonal differences).                  |
# | **1** | **Low Severity**                                      | Light surface burn: understory burned, but canopy intact; minimal mortality.                               |
# | **2** | **Moderate Severity**                                 | Some canopy scorching and mortality; noticeable vegetation change; soil impacts minor-to-moderate.         |
# | **3** | **High Severity**                                     | Very strong vegetation change; canopy consumption or full mortality; strong soil and ground cover impacts. |
# | **4** | **Increased Greenness / Regrowth (sometimes 4 = NA)** | Area shows greater greenness post-fire than pre-fire; often mesic sites or recovering burned areas.        |
# | **5** | **Non-Processing Area Mask (NPAM)**                   | Clouds, shadows, snow, water, or other regions where severity couldn't be calculated.                      |
# | **6** | **Missing / No Data**                                 | Outside MTBS processing region or inaccessible imagery.                                                    |

recode_severity <- function(r) {
  # reclass matrix: from, to, becomes
  m <- matrix(c(
    0, 0, 0,   # unburned/regrowth -> 0
    1, 1, 1,   # low -> 1
    2, 2, 1,   # moderate -> 1
    3, 3, 2,   # high -> 2
    4, 6, NA   # classes 4–6 → NA
  ), ncol = 3, byrow = TRUE)

  classify(r, m, include.lowest = TRUE)
}

# Apply to each raster
recode_rasters <- lapply(rasters, recode_severity)

ex <- rasters[[1]]
ex1 <- ex[[1]]
ex1
png(here(plotpath, "burn_severity_example.png"), width = 6, height = 6, units = "in", res = 300)
plot(ex1)
dev.off()
