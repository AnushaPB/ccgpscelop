get_feems <- function(){
  # Code from this issue: https://github.com/NovembreLab/feems/issues/34
  # Read in output from FEEMS
  feems_nodes <- read_csv(here("analysis", "feems", 'feems_nodes.csv'), col_names='node_id', col_types='i')
  feems_node_pos <- read_csv(here("analysis", "feems", 'feems_node_pos.csv'), col_names=c('lon', 'lat'), col_types='dd')
  feems_edges <- read_csv(here("analysis", "feems", 'feems_edges.csv'), col_names=c('n1', 'n2'), col_types='ii')
  feems_w <- read_csv(here("analysis", "feems", 'feems_w.csv'), col_names='w', col_types='d')

  # Add lon/lat to nodes and weights to edges
  nodes <- add_column(feems_nodes, feems_node_pos)
  edges <- add_column(feems_edges, feems_w) 

  # Join nodes and weights to edges
  edges_n1 <- left_join(edges, nodes, by=c('n1' = 'node_id'))
  edges_n2 <- left_join(edges_n1, nodes, by=c('n2' = 'node_id'))
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

  # Crop to CA
  if (crop) {
    ca <- get_ca()
    feems <- st_intersection(feems, ca)
  }

  # Add nodes
  if (nodes){
    nodes_sf <- 
      nodes %>%
      st_as_sf(coords = c("lon", "lat"), crs = 4326)
    feems <- list(edges = feems, nodes = nodes_sf)
  }
  
  return(feems)
}