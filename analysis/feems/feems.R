get_feems <- function(crop = TRUE, nodes = FALSE){
  
  # Read in edges
  feems <- st_read(here("analysis", "feems", "feems_edges.shp"))

  # Assign CRS
  st_crs(feems) <- st_crs(4326)
  
  # Log transform like FEEMS does to get comparison to mean
  feems$log_weight <- log10(feems$weight) - mean(log10(feems$weight), na.rm = TRUE)
  
  # Crop to CA
  if (crop) {
    ca <- get_ca()
    feems <- st_intersection(feems, ca)
  }

  # Add nodes
  if (nodes){
    nodes <- st_read(here("analysis", "feems", "feems_nodes.shp"))
    st_crs(nodes) <- st_crs(4326)
    feems <- list(edges = feems, nodes = nodes)
  }

  return(feems)
}
