
# functions to return and calculate bias layers
get_bias <- function(spp, rank, n = 10000, lyr = NULL, limit = 1000, cache = TRUE, dir = "outputs"){
  bb <- rgbif::name_backbone(spp)
  
  if(rank == "phylum") {
    key <- bb$phylumKey
    bias_path <- here(dir, paste0(key, "_", rank, "_", n, "_bias.csv"))
  }
  
  if(rank == "class") {
    key <- bb$classKey
    bias_path <- here(dir, paste0(key, "_", rank, "_", n, "_bias.csv"))
  }

  if (file.exists(bias_path) & cache) {
    bias_bkg <- read.csv(bias_path)
  } else {
    if (is.null(lyr)) "lyr must be provided"
    bias_bkg <- make_biasbkg(key, rank, lyr, limit, n)
  }
  
  return(bias_bkg)
}

get_lyr <- function(){
  file_path <- here("data", "env", "paleoclim", "CHELSA_cur_V1_2B_r2_5m.zip")

  # Load and process the data using purrr
  ca <- get_ca()
  data <- rpaleoclim::load_paleoclim(file_path)
  data <- crop(data, ca)
  data <- mask(data, ca)

  lyr <- data[[1]]*0
  return(lyr)
}

make_biasbkg <- function(key, rank, lyr, limit = 1000, n = 10000){
  # get key for class or phylum
  if (is.character(key) & rank == "phylum") key <- rgbif::name_backbone(key)$phylumKey
  if (is.character(key) & rank == "class") key <- rgbif::name_backbone(key)$classKey
  
  # get GBIF points
  bias_dat <- get_biasgbif(key, rank, limit)

  # make into raster
  bias_lyr <- make_biaslyr(bias_dat, lyr)

  # save output
  terra::writeRaster(bias_lyr, here("outputs", paste0(key, "_", rank, "_", n, "_bias.tif")), overwrite = TRUE)
  
  # get coords for background points
  bias_bkg <- sample_bias(bias_lyr, n)
  write.csv(bias_bkg, here("outputs", paste0(key, "_", rank, "_", n, "_bias.csv")), row.names = FALSE)
  
  message(paste(rank, key, "complete"))
  
  return(bias_bkg)
}


get_biasgbif <- function(key, rank, limit = 1000){
  
  if (rank == "phylum"){
    gbif <- rgbif::occ_search(phylumKey = key, hasCoordinate = T, stateProvince = "California", limit = limit)
  }
  if (rank == "class"){
    gbif <- rgbif::occ_search(classKey = key, hasCoordinate = T, stateProvince = "California", limit = limit)
  }
  
  bias_dat <- gbif$data
  return(bias_dat)
}

make_biaslyr <- function(bias_dat, lyr, agg = NULL){
  # Create sf coordinates
  bias_dat <- st_as_sf(bias_dat, coords = c("decimalLongitude", "decimalLatitude"))
  
  # aggregate raster layer
  if (!is.null(agg)) lyr <- terra::aggregate(lyr, agg)
  
  # Create raster layer of counts for each cell
  taxon_bias = terra::rasterize(bias_dat, lyr, fun = "count", background = NA)
  
  # Assign NA values of bias raster to the min
  bias_mean = taxon_bias
  bias_mean[is.na(bias_mean)] = terra::global(bias_mean, "min", na.rm = TRUE)
  bias_lyr = mask(bias_mean, lyr)
  
  return(bias_lyr)
  
}

sample_bias <- function(bias_lyr, n = 10000){
  # Convert from count to probability
  bias_lyr_probs <- bias_lyr / sum(values(bias_lyr), na.rm = TRUE)

  # Replace NAs in probability vector with 0s (don't want to sample from NA areas)
  bias_lyr_probs[is.na(bias_lyr_probs)] <- 0

  # Generate a random sample of indices based on bias/probability
  sample_indices <- sample(x = 1:ncell(bias_lyr), size = n, replace = TRUE, prob = values(bias_lyr_probs))

  # Extract the sampled values (or coordinates) from the raster
  # This gives you the cell values; for coordinates, use xyFromCell
  bias_bkg <- terra::extract(bias_lyr, sample_indices, xy = TRUE)[,c("x", "y")]

  return(bias_bkg)
}

allspp_bias <- function(spps, lyr, rank, limit = 1000, n = 10000){
  bb <- purrr::map_dfr(spps, name_backbone)
  
  phylkey <- unique(bb$phylumKey)
  clsskey <- unique(bb$classKey)
  
  lsp <- sapply(phylkey, make_biasbkg, rank = "phylum", lyr = lyr, limit = limit, n = n)
  lsc <- sapply(clsskey, make_biasbkg, rank = "class", lyr = lyr, limit = limit, n = n)
  
  return(list(phylbias = lsp, clssbias = lsc))
}
