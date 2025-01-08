get_ca <- function() {
  # Load the U.S. state boundaries data
  states <- tigris::states(cb = TRUE)
  # Extract the boundary of California (CA)
  ca <- states[states$STUSPS == "CA", "STUSPS"]
  ca <- sf::st_transform(ca, sf::st_crs(4326))
  return(ca)
}

get_coords <- function(sf = FALSE) {
  # sample coords
  coords <- read_table(here("data", "ccgp_data", "58-Sceloporus.coords.txt"), col_names = FALSE)
  colnames(coords) <- c("SampleID", "x", "y")
  if (sf) {
    coords <- st_as_sf(coords, coords = c("x", "y"), crs = 4326)
  }
  return(coords)
}

get_dem <- function(r = FALSE) {
  coords <- get_coords(sf = TRUE)
  ca <- get_ca()
  dem <- elevatr::get_elev_raster(coords, z = 5)
  dem <- mask(dem, ca)
  dem <- crop(dem, ca)
  if (r) return(dem)
  coords_dem <- st_transform(coords, st_crs(dem))
  dem_df <- terra::extract(dem, coords_dem, ID = FALSE)
  dem_df <- data.frame(SampleID = coords$SampleID, elev = dem_df)
  return(dem_df)
}

get_env <- function() {
  read_csv(here("data", "env", "envdata.csv"))
}

get_biokey <- function(){
  biokey <- 
    data.frame(
      bio = paste0("bio_", 1:20),
      description = c(
        "mean annual temperature",
        "min temperature of coldest month",
        "max temperature of warmest month",
        "temperature annual range",
        "mean diurnal range",
        "isothermality",
        "temperature seasonality",
        "max temperature of warmest period",
        "min temperature of coldest period",
        "temperature annual range",
        "mean temperature of wettest quarter",
        "mean temperature of driest quarter",
        "mean temperature of warmest quarter",
        "mean temperature of coldest quarter",
        "annual precipitation",
        "precipitation of wettest month",
        "precipitation of driest month",
        "precipitation seasonality",
        "precipitation of wettest quarter",
        "precipitation of driest quarter"
      )
    ) %>%
    mutate(
      vartype = 
        case_when(
          grepl("temperature", description) ~ "temperature", 
          grepl("precipitation", description) ~ "precipitation",
          grepl("isothermality", description) ~ "temperature",
          grepl("mean diurnal range", description) ~ "temperature",
          TRUE ~ description
          )
        )
    
    return(biokey)
  }

get_range <- function(){
  sf::st_read(here("data", "rWFLIx_CONUS_HabMap_2001v1")) %>% 
    st_transform(3310)
}
