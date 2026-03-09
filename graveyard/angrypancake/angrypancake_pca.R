library(tidyverse)
library(algatr)
library(here)
library(raster)
library(terra)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

datapath <- here("analysis", "angrypancake", "data")
plotpath <- here("analysis", "angrypancake", "plots")

# IDEA ADD REDUNDANT ALLELES?
coords <- 
  read_table(here(datapath, "58-Sceloporus.coords.txt"), col_names = c("SampleID", "x", "y")) %>%
  st_as_sf(coords = c("x", "y"), crs = 4326) %>%
  st_transform(3310)

na <- ne_countries(scale = "medium", returnclass = "sf") 
na_coords <- coords %>% st_buffer(500000) %>% st_transform(st_crs(na))
na <- st_crop(st_make_valid(na), st_bbox(na_coords)) %>% st_transform(3310)

gen <- read_table(here(datapath, "nonsyn.raw"))
dos <- 
  gen %>% 
  dplyr::select(SampleID = IID, starts_with("chr")) 

coords <- coords %>% filter(SampleID %in% dos$SampleID)

# Confirm order of samples is the same
stopifnot(all(coords$SampleID == dos$SampleID))

dos <- dos %>% dplyr::select(-SampleID)
allele <- dos
allele[allele > 0] <- 1

# current_env <- rast(here(datapath, "CHELSA_bio1_1981-2010_V.2.1.tif"))
# names(current_env) <- "bio1"
# future_env <- rast(here(datapath, "CHELSA_bio1_2071-2100_gfdl-esm4_ssp585_V.2.1.tif"))
# # names(future_env) <- "bio1"
# latlon <- na %>% st_transform(4326)
# current <- crop(current_env, latlon)
# future <- crop(future_env, latlon)
# current <- project(current, "EPSG:3310")
# future <- project(future, "EPSG:3310")
# TEMPORARY
current <- aggregate(current, 2)
future <- aggregate(future, 2)
current <- mask(current, na)
future <- mask(future, na)
# BIOCLIM DATA
global_bioclim <- rast(list.files(here("data", "env", "chelsa"), full.names = TRUE))
ca_bioclim <- crop(global_bioclim, st_transform(na, 4326))
current <- ca_bioclim %>% project("epsg:3310")
# Rename bioclim
names(current) <- gsub("_1981-2010_V.2.1", "", gsub("CHELSA_", "", names(current)))
# Read in future environmental layers and project to 3310
future_paths <- list.files(here("data", "env", "chelsa_future"), full.names = TRUE)
future_paths <- future_paths[grepl("bio", future_paths)]  # only bioclim
global_future <- rast(future_paths)
ca_future <- crop(global_future, st_transform(na, 4326))
future <- ca_future %>% project("epsg:3310")
# Rename future layers
names(future) <- gsub("_2071-2100_gfdl-esm4_ssp585_V.2.1", "", gsub("CHELSA_", "", names(future)))
# SUBSET ENVSTACK TO MATCH FUTURE STACK
# Some bioclim variables are missing in the future dataset, so we subset envstack to only those variables that are present in the future dataset
current <- current[[names(future)]]
all(names(current) ==  names(future))  # should be 0
# PCA FUNCTION
run_pca <- function(cur, fut){
  # Run a RasterPCA on model envstack
  # 1) Compute scaling from CURRENT layers
  mu  <- global(cur, "mean", na.rm = TRUE)$mean
  sig <- global(fut, "sd",   na.rm = TRUE)$sd

  # 2) Get values (matrix: ncell x nlayers) and fit PCA using current's mu/sig
  #    prcomp will do the scaling internally when we pass center/scale.
  # center: a logical value indicating whether the variables should be
  #         shifted to be zero centered. Alternately, a vector of length
  #         equal the number of columns of ‘x’ can be supplied.  The
  #         value is passed to ‘scale’.

  # scale.: a logical value indicating whether the variables should be
  #         scaled to have unit variance before the analysis takes place.
  #         The default is ‘FALSE’ for consistency with S, but in general
  #         scaling is advisable.  Alternatively, a vector of length
  #         equal the number of columns of ‘x’ can be supplied.  The
  #         value is passed to ‘scale’.
  X <- values(cur, mat = TRUE)
  cc <- complete.cases(X)
  pca <- prcomp(X[cc, , drop = FALSE], center = mu, scale. = sig, retx = FALSE)

  # 3) Project the first 3 PCs for current & future with the SAME PCA object
  model_current <- predict(cur, pca, index = 1:3)
  model_future  <- predict(fut, pca, index = 1:3)

  return(list(model_current = model_current, model_future = model_future))
}
pca_result <- run_pca(current, future)
current_pc <- pca_result$model_current
future_pc <- pca_result$model_future

# Extract environmental vars
env <- raster::extract(current_pc, coords, ID =FALSE)

# Standardize environmental variables and make into dataframe
env <- scale(env, center = TRUE, scale = TRUE)
env <- data.frame(env)

# plot(current)
# points(coords)
tp <- function(x) {
  png(here("TEMP.png"))
  plot(x)
  dev.off()
}

# GO OVER MAP1/2 thing (maybe change ot A/B for clarity)
predict_allele <- function(snp, coords, current, future, maxdist = 50000, admix_thresh = 0) {
  current_vals <- extract(current, coords, ID = FALSE)
  future_vals <- extract(future, coords, ID = FALSE)

  # Model current allele presence/absence
  current_df <- data.frame(allele_ref = snp %in% c(0, 1), allele_alt = snp %in% c(1, 2), current_vals)
  mod_ref <- glm(allele_ref ~ PC1 + PC2 + PC3, data = current_df, family = binomial)
  mod_alt <- glm(allele_alt ~ PC1 + PC2 + PC3, data = current_df, family = binomial)

  # Predict current distribution
  current_map_ref <- terra::predict(current, mod_ref, type = "response")
  current_map_ref <- ifel(current_map_ref < 0, 0, current_map_ref)
  current_map_ref <- ifel(current_map_ref > 1, 1, current_map_ref)
  current_map_bin_ref <- current_map_ref > 0.5
  current_map_alt <- terra::predict(current, mod_alt, type = "response")
  current_map_alt <- ifel(current_map_alt < 0, 0, current_map_alt)
  current_map_alt <- ifel(current_map_alt > 1, 1, current_map_alt)
  current_map_bin_alt <- current_map_alt > 0.5

  # Predict which allele needs to be present in the future
  future_map_ref <- terra::predict(future, mod_ref, type = "response")
  future_map_ref <- ifel(future_map_ref < 0, 0, future_map_ref)
  future_map_ref <- ifel(future_map_ref > 1, 1, future_map_ref)
  future_map_bin_ref <- future_map_ref > 0.5
  future_map_alt <- terra::predict(future, mod_alt, type = "response")
  future_map_alt <- ifel(future_map_alt < 0, 0, future_map_alt)
  future_map_alt <- ifel(future_map_alt > 1, 1, future_map_alt)
  future_map_bin_alt <- future_map_alt > 0.5

  # Combine according to conditions:
  # (ref, !alt) → 0
  # (!ref, alt) → 1
  # (ref, alt) → 0.5
  current_map_bin <- cover(
    ifel(current_map_bin_ref & !current_map_bin_alt, 0,    NA),
    ifel(!current_map_bin_ref &  current_map_bin_alt, 1,    NA)
  )
  # Deal with hets: if both are TRUE, set to 0.5, then see if one is clearly more likely than the other (by > 0.1), if so set to that allele, otherwise keep as 0.5
  current_map_bin[current_map_bin_ref & current_map_bin_alt] <- 0.5
  current_map_bin[current_map_bin == 0.5 & (current_map_ref - current_map_alt) > admix_thresh] <- 0
  current_map_bin[current_map_bin == 0.5 & (current_map_alt - current_map_ref) > admix_thresh] <- 1

  future_map_bin <- cover(
    ifel(future_map_bin_ref & !future_map_bin_alt, 0,    NA),
    ifel(!future_map_bin_ref &  future_map_bin_alt, 1,    NA)
  )
  # Deal with hets: if both are TRUE, set to 0.5, then see if one is clearly more likely than the other (by > 0.1), if so set to that allele, otherwise keep as 0.5
  future_map_bin[future_map_bin_ref & future_map_bin_alt] <- 0.5
  future_map_bin[future_map_bin == 0.5 & (future_map_ref - future_map_alt) > admix_thresh] <- 0
  future_map_bin[future_map_bin == 0.5 & (future_map_alt - future_map_ref) > admix_thresh] <- 1
 
  #plot(future_map_bin
  # Get difference (and treat 0.5 as no difference)
  diff_map <- ifel(
    # If the same, no difference
    future_map_bin == current_map_bin, 0,
    # If either are 0.5, no difference, if not → difference
    ifel(future_map_bin == 0.5 | current_map_bin == 0.5, 0, 1)
  )
  # Get raw difference
  raw_diff_map <- future_map_bin != current_map_bin
  # Get difference between raw_diff_map and diff_map to quantify how much heterozygosity is the saving difference
  # diff_map (A) will always be < raw_diff_map (B)
  # 0 = no change, 1 = change
  # Possible values: B = 0, A = 0 → 0 (no change, no difference)
  #                  B = 1, A = 0 → 1 (change, and difference due to het)
  #                  B = 1, A = 1 → 0 (change, and no difference due to het)
  het_diff_map <- raw_diff_map - diff_map
  #plot(diff_map)

  # ADAPTING
  # Targets for distance(): set target cells to 1 and everything else to NA
  cur_ref_targets <- ifel(current_map_bin_ref == 1, 1, NA)   # where allele ref is present now
  cur_alt_targets <- ifel(current_map_bin_alt == 1, 1, NA)   # where allele alt is present now

  # Distance (map units) to nearest needed-allele cell in the CURRENT map
  # MASK to current map so distances are only calculated for land cells
  d_current_ref <- mask(distance(cur_ref_targets), current[[1]])
  d_current_alt <- mask(distance(cur_alt_targets), current[[1]])

  # Combine
  adapt_distance <- 
    # If the allele persists at the same cell or is 0.5 → distance = 0
    ifel(current_map_bin == future_map_bin | current_map_bin == 0.5 | future_map_bin == 0.5, 0,
      # Otherwise → distance to nearest CURRENT occurrence of that allele
      # In other words: if future needs ref (0), then get the distance in the current to 0. The else just means that if future needs alt (1), get the distance in the current to 1
      ifel(future_map_bin == 0, d_current_ref, d_current_alt))
  adapt_met <- adapt_distance <= maxdist

  # MOVING
  # Future targets: where each allele will exist in the FUTURE
  fut_ref_targets <- ifel(future_map_bin_ref == 1, 1, NA)
  fut_alt_targets <- ifel(future_map_bin_alt == 1, 1, NA)

  # Distances (map units; meters if CRS is projected) to nearest FUTURE cell with that allele
  d_future_ref <- mask(distance(fut_ref_targets), future[[1]])
  d_future_alt <- mask(distance(fut_alt_targets), future[[1]])

  # Movement needed to keep the same allele (track forward):
  # - 0 if current allele is 0.5
  # - 0 if the allele persists in place (because that cell is also a future target)
  # - else distance to nearest FUTURE location holding that allele
  # CHECK THIS
  move_distance <- 
    # If the allele persists at the same cell or is 0.5 → distance = 0
    ifel(current_map_bin == future_map_bin | current_map_bin == 0.5 | future_map_bin == 0.5, 0,
        # Otherwise → distance to nearest FUTURE occurrence of that allele
        # In other words: if current is ref (0), then get the distance in the future to 0. The else just means that if current is alt (1), get the distance in the future to 1
        ifel(current_map_bin == 0, d_future_ref, d_future_alt))
  move_met <- move_distance <= maxdist
  
  result <- c(adapt_distance, adapt_met, move_distance, move_met, diff_map, raw_diff_map, het_diff_map, current_map_bin, future_map_bin)
  names(result) <- c("adapt_distance", "adapt_bin", "move_distance", "move_bin", "diff", "raw_diff", "het_diff", "current_bin", "future_bin")
  result <- as.list(result)
  return(result)
}

snp_names <- colnames(dos)
names(snp_names) <- snp_names
current_agg <- aggregate(current_pc, 5)
future_agg <- aggregate(future_pc, 5)
set.seed(9204)
availability <- 
  map(sample(snp_names, 50), ~{
    #.x = "chr1_192013_T_A_T"
    snp <- dos %>% pull(.x)
    predict_allele(snp, coords, current_agg, future_agg, maxdist = 10000, admix_thresh = 0.1)
  }, .progress = TRUE) 

ca <- 
  rnaturalearth::ne_states(country = "United States of America", returnclass = "sf") %>%
  filter(name == "California") %>%
  st_transform(3310)

current_ca <- mask(current_agg, ca)
avail_crop <- map(availability, ~as.list(crop(mask(rast(.x), current_ca), ca)), .progress = TRUE)

avail_stack <- list_transpose(avail_crop)
names(avail_stack) <- names(rast(availability[[1]]))

avail_stack <- map(avail_stack, rast)
mean_stack <- map(avail_stack, ~mean(.x), na.rm = TRUE) %>% rast()
dist_stack <- mean_stack[[c("adapt_distance", "move_distance")]]
dist_stack[dist_stack == 0] <- NA

gg_df <- 
  extract(mean_stack, coords, ID = FALSE) %>% 
  as.data.frame() %>%
  bind_cols(coords) %>%
  st_as_sf()

plot_map <- function(var){
  ggplot() +
    geom_sf(data = ca, fill = "lightgrey", color = "darkgrey") +
    geom_sf(data = gg_df, aes(col = .data[[var]])) +
    scale_color_viridis_c(option = "plasma", direction = 1, na.value = "transparent") +
    coord_sf() +
    theme_void() +
    labs(col = var)
}

pdf(here(plotpath, "angrypancake_maps.pdf"), width = 12, height = 10)
plot(mean_stack, col = viridis::plasma(100, direction = 1))
plot(dist_stack, col = viridis::plasma(100, direction = 1))
plot(avail_stack$diff, col = viridis::plasma(100, direction = 1))
plot(avail_stack$raw_diff, col = viridis::plasma(100, direction = 1))
plts <- map(names(mean_stack), ~(plot_map(.x)))
cowplot::plot_grid(plotlist = plts, ncol = 3)
dev.off()

# Plot example
example1 <- avail_crop[[1]] %>% rast()
snp1 <- dos %>% pull(colnames(dos)[1]) 

example2 <- avail_crop[[2]] %>% rast()
snp2 <- dos %>% pull(colnames(dos)[2]) 

example3 <- avail_crop[[3]] %>% rast()
snp3 <- dos %>% pull(colnames(dos)[3])

pdf(here(plotpath, "angrypancake_example.pdf"), width = 10, height = 10)
plot(example1)
plot(example2)
plot(example3)
dev.off()

# ADD IN HOW FAR THEY WOULD HVAE TO MOVE 

# NEED TO MAKE MAPS OF ACTUAL ALLELE FREQUENCIES NOT JUST PREDICTED ALLELE FREQUENCIES (SINCE DOESNT REALLY WORK WELL FOR HETEROZYGOTES)
# MAYBE TRYE MORE ENV VALUES

# ACTUALLY, CANT USE KRIGIN SINCE IT IS SO SAMPLING DPEENDENT \\


# ADD SOEEMTHING TO CALCULATE WHERE THE SOURCE OF ALLELES ARE FOR ADAPT AND WHERE THE MOVE DESTINATIONS ARE FOR MOVE