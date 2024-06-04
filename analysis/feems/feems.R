get_feems <- function(crop = TRUE, nodes = TRUE){
  # Code from this issue: https://github.com/NovembreLab/feems/issues/34
  # Read in output from FEEMS
  output_dir <- here("analysis", "feems", "outputs")
  feems_nodes <- read_csv(here(output_dir, 'feems_nodes.csv'), col_names='node_id', col_types='i')
  feems_node_pos <- read_csv(here(output_dir, 'feems_node_pos.csv'), col_names=c('lon', 'lat', "nsamp"), col_types='dd')
  feems_edges <- read_csv(here(output_dir, 'feems_edges.csv'), col_names=c('n1', 'n2'), col_types='ii')
  feems_w <- read_csv(here(output_dir, 'feems_w.csv'), col_names='w', col_types='d')

  # Add lon/lat to nodes and weights to edges
  feems_nodes <- add_column(feems_nodes, feems_node_pos)
  feems_edges <- add_column(feems_edges, feems_w) 

  # Join nodes and weights to edges
  edges_n1 <- left_join(feems_edges, feems_nodes, by=c('n1' = 'node_id'))
  edges_n2 <- left_join(edges_n1, feems_nodes, by=c('n2' = 'node_id'))
  feems_edges <- select(edges_n2, lon1=lon.x, lat1=lat.x, lon2=lon.y, lat2=lat.y, w)

  # Convert edges to sf object
  feems <- 
    feems_edges %>%
    mutate(
      geometry = pmap(list(lon1, lat1, lon2, lat2), \(lon1, lat1, lon2, lat2) {
        st_linestring(matrix(c(lon1, lat1, lon2, lat2), ncol = 2, byrow = TRUE))
      })
    ) %>%
    st_as_sf() %>%
    select(-lon1, -lat1, -lon2, -lat2) %>%
    st_set_crs(4326)

  # Make transformed weight
  feems$w_trans <- log10(feems$w) - mean(log10(feems$w), na.rm = TRUE)
  
  # Crop to CA
  if (crop) {
    ca <- get_ca()
    feems <- st_intersection(feems, ca)
  }

  # Add nodes
  if (nodes){
    nodes_sf <- 
      feems_nodes %>%
      st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
      mutate(nsamp = case_when(nsamp == 0 ~ NA, TRUE ~ nsamp)) 
    feems <- list(edges = feems, nodes = nodes_sf)
  }
  
  return(feems)
}
