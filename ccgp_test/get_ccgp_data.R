
get_ccgp_het <- function(){
  het_files <- list.files("~/ccgp/het", pattern = ".het", full.names = TRUE)
  nsites <- read_csv("~/ccgp/het/bed_nsites.csv", col_names = c("project", "nsites"))
  het <- 
    purrr::map(het_files, ~{
      # Read in data
      het <- readr::read_table(.x)

      # Calculate number of homozygotes 
      het <- het %>% mutate(Ho = (`N(NM)` - `O(HOM)`))

      if ("IID" %in% names(het)) het <- rename(het, SampleID = IID)
      if ("INDV" %in% names(het)) het <- rename(het, SampleID = IID)

      # Add project ID
      het <- mutate(het, project = stringr::str_extract(.x, "\\d+-\\w+"))

      return(het)

    }) %>% 
    bind_rows()

  # Add number of sites
  het_df <- 
    left_join(het, nsites, by = "project") %>%
    mutate(Ho = Ho/nsites)

  return(het_df)
}


get_ccgp_coords <- function(sf = TRUE){
  coords_files <- list.files("~/ccgp/coords", pattern = ".coords", full.names = TRUE)
  coords <- 
    purrr::map(coords_files, ~{
      # Read in data
      coords <- readr::read_table(.x, col_names = FALSE)

      # Name the columns
      coords <- rename(coords, SampleID = X1, x = X2, y = X3)

      # Add project ID
      coords <- mutate(coords, project = stringr::str_extract(.x, "\\d+-\\w+"))

      return(coords)

    }) %>% 
    bind_rows()

  if (sf) {
    coords <- st_as_sf(coords, coords = c("x", "y"), crs = 4326) %>% st_transform(3310)
  }
}
