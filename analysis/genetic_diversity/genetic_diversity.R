
get_het <- function(){
  # callable sites for denominator
  callable <- read.csv(here("data_processing", "callable_counts.csv"))
  callable_sites <- callable %>% pull(callable_sites)

  het <- format_het(here("data", "ccgp_data", "58-Sceloporus.het"), callable_sites = callable_sites)
  
  if ("IID" %in% names(het)) het$SampleID <- het$IID
  if ("INDV" %in% names(het)) het$SampleID <- het$INDV

  return(het)
}

format_het <- function(path, callable_sites){
  # Load the data in R
  het_data <- read_table(path)

  # Calculate the average heterozygosity per individual
  if ("N(NM)" %in% names(het_data)) het_data$Ho <- (het_data$`N(NM)` - het_data$`O(HOM)`)/callable_sites
  if ("N_SITES" %in% names(het_data)) het_data$Ho <- (het_data$N_SITES - het_data$`O(HOM)`)/callable_sites
  return(het_data) 
}

get_pops <- function(){
  map(paste0("pop", 1:5), ~{
    read.table(here("analysis", "admixture", "outputs", paste0("k5_", .x, ".txt")), header = FALSE) %>%
    mutate(pop = .x) %>%
    rename(SampleID = V1)
  }) %>% bind_rows()
}

get_roh <- function(){
  coords <- get_coords(sf = TRUE)
  pops <- get_pops()
  coords <- left_join(coords, pops, by = "SampleID")

  roh <-
    map(paste0("pop", 1:5), ~{
      poproh <- 
        read_table(here("analysis", "genetic_diversity", "outputs", paste0("K5_", .x, ".froh")), col_names = FALSE) %>%
        rename(SampleID = X1, froh = X2) %>%
        right_join(filter(coords, pop == .x)) %>%
        # missing froh values indicates no roh greater than the minimum size were found
        # remains as NA if pop was not found (i.e., individual was not included in analysis)
        mutate(
          froh0 = case_when(is.na(froh) & !is.na(pop) ~ 0, .default = froh),
          pop = .x
        )
    }) %>% bind_rows()

  return(roh)
}


geosummarize <- function(coords, stat = "Ho", res = 50000){
  # transform coords
  coords <- coords %>% st_transform(3310)

  # Load or create a raster
  lyr <- wingen::coords_to_raster(coords, res = res)

  # Get raster coords
  raster_coords <- terra::extract(lyr, coords, xy = TRUE)

  # combine with coords
  final_coords <- 
    bind_cols(raster_coords, coords) %>%
    group_by(x, y) %>%
    summarize_at({{stat}}, mean, na.rm = TRUE) %>%
    drop_na(x, y) %>%
    st_as_sf(coords = c("x", "y"), crs = st_crs(coords)) %>%
    st_transform(4326)

}


