
get_ccgp_het <- function(){
  het_files <- list.files("~/ccgp/het", pattern = ".het", full.names = TRUE)
  het <- 
    purrr::map(het_files, ~{
      # Read in data
      het <- readr::read_table(.x)

      # Calculate heterozygosity 
      # TODO: CHANGE DENOMINATOR TO BED FILE LENGTH
      het <- het %>% mutate(Ho = (`N(NM)` - `O(HOM)`)/`N(NM)`)

      if ("IID" %in% names(het)) het <- rename(het, SampleID = IID)
      if ("INDV" %in% names(het)) het <- rename(het, SampleID = IID)

      # Add project ID
      het <- mutate(het, project = stringr::str_extract(.x, "\\d+-\\w+"))

      return(het)

    }) %>% 
    bind_rows()
}


get_ccgp_coords <- function(){
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
}

get_ccgp_taxa <- function(){
  # Get coordinates
  df <- get_ccgp_coords()

  # Get list of species names
  taxa_df <- 
    data.frame(project = df$project) %>%
    mutate(genus = gsub("[^A-Za-z]", "", project)) %>%
    distinct() %>%
    drop_na(project)

  # Get taxonomic information
  taxize::taxize_options(ncbi_sleep = 0.4)
  taxa_info <- taxize::classification(taxa_df$genus, db = "ncbi")

  # Convert to df
  taxa_ls <- map(1:length(taxa_info), ~taxa_info[[.x]])
  names(taxa_ls) <- taxa_df$project

  # Get the values with NA
  taxa_null <- map(taxa_ls, ~if(length(.x) == 1) return(NULL) else return(.x))
  which(is.null(taxa_null))

  taxa_df <- 
    taxa_null %>%
    compact() %>%
    bind_rows(.id = "project") %>%
    filter(rank %in% c("kingdom", "phylum", "class", "order", "family", "genus", "species"))

  return(taxa_df)

}
