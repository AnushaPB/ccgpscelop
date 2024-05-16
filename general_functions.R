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
  coords <- read_table(here("data", "raw_data", "58-Sceloporus.coords.txt"), col_names = FALSE)
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
  dem <- as.data.frame(dem, xy = TRUE)
  colnames(dem) <- c("x", "y", "elev")
  return(dem)
}

get_env <- function() {
  read_csv(here("data", "env", "envdata.csv"))
}
