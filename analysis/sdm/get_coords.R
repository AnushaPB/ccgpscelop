# Functions to get and clean GBIF and CCGP coords

#' CCGP / GBIF coordinate cleaner
#' @return cleaned coordinates to use as a 'create_sdm' argument
#' 
# TODO: CHECK AND CLEAN THIS FUNCTION
get_occ_coords <- function(spp, data_source = "both", occ_limit = 10000, envlayers = NULL, cache = TRUE, stateProvince = "California", grid_sample = TRUE, ncores = NULL){
  
  if (data_source == "ccgp" | data_source == "both") ccgp <- get_coords() %>% mutate(decimalLongitude = x, decimalLatitude = y)
  
  #if not combining with gbif data, clean data is only ccgp data
  if (data_source == "ccgp") dat_cl <- ccgp
 
  if (data_source == "both" | data_source == "gbif") dat_cl <- get_gbif(spp, occ_limit = occ_limit, stateProvince = stateProvince, ncores = ncores, cache = cache)
    
  if (data_source == "both") dat_cl <- bind_rows(dat_cl, ccgp)  

  coords <- dat_cl %>% dplyr::select(decimalLongitude, decimalLatitude) 
  
  if (grid_sample & !is.null(envlayers)) coords <- dismo::gridSample(coords, r =envlayers, n = 1)
  
  return(coords)
}


# Helper function to get CCGP coordinates
get_coords_ccgp <- function(spp, ID = NULL, keepID = FALSE, AddSppCol = FALSE){
  # get full WGS data
  gsd <- read_csv(here("data", "WGS_METADATA_DB.csv"), show_col_types = FALSE)
  
  # get species data
  df <- gsd[grepl(spp, gsd$`*organism`), ]
  
  # get coords
  coords <- df %>% 
    dplyr::select(`*sample_name`, long, lat) %>% 
    mutate_at(c("long", "lat"), as.numeric) %>%
    rename("ID" = `*sample_name`, "decimalLongitude" = "long", "decimalLatitude" = "lat") 
  
  # subset IDs
  if(!is.null(ID)){coords <- coords[coords$ID %in% ID,]}
  
  # remove ID 
  if(!keepID){coords <- coords %>% dplyr::select(decimalLongitude, decimalLatitude)}
  
  # add spp
  if(AddSppCol){coords$species <- spp}
  
  return(coords)
}

get_gbif <- function(spp, occ_limit = 10000, stateProvince = "California", ncores = NULL, cache = TRUE){
  
  spp <- gsub(" ", "_", spp)
  if (is.null(stateProvince)) {
    path <- here("data", "gbif", paste0(spp, "_gbif_global.csv"))
  } else if (stateProvince == "California") path <- here("data", "gbif", paste0(spp, "_gbif_CA.csv"))
  spp <- gsub("_", " ", spp)
  
  if (cache & file.exists(path)) {
    message("Using existing GBIF file")
    return(read.csv(path))
  }
  
  if (is.null(ncores)) {
    dat <- 
      purrr::map(seq(0, occ_limit, 1000), spp_search, spp = spp, stateProvince = stateProvince, .progress = TRUE) %>%
      bind_rows() %>%
      distinct()
  } else {
    future::plan(future::multisession, workers = ncores)
    
    dat <- 
      furrr::future_map(seq(0, occ_limit, 1000), spp_search, spp = spp, stateProvince = stateProvince, 
                        .options = furrr::furrr_options(seed = 20), .progress = TRUE) %>%
      bind_rows() %>%
      distinct()

    future::plan("sequential")
  }

  if (nrow(dat) == 0) {warning("no occurrence records found, returning NULL"); return(NULL)}
  
  dat_cl <- cleaner(dat)
  
  if (cache) write.csv(dat_cl, path, row.names = FALSE)

  return(dat_cl)
}

cleaner <- function(dat){
  if (nrow(dat) == 0) {warning("data has 0 rows, returning NULL"); return(NULL)}
  
  dat$countryCode <- countrycode(dat$countryCode, origin =  'iso2c', destination = 'iso3c')
  
  dat <- data.frame(dat)
  
  flags <- clean_coordinates(x = dat,
                             lon = "decimalLongitude",
                             lat = "decimalLatitude",
                             countries = "countryCode",
                             species = "species",
                             tests = c("capitals", "centroids", "equal","gbif", "institutions",
                                       "zeros", "countries")) # most test are on by default
  
  #Exclude problematic records
  dat_cl <- dat[flags$.summary,]
  
  # https://ropensci.github.io/CoordinateCleaner/articles/Cleaning_GBIF_data_with_CoordinateCleaner.html
  # cf_age is used for fossil cleaning and removing temporal outliers
  flags <- cf_age(x = dat_cl,
                  lon = "decimalLongitude",
                  lat = "decimalLatitude",
                  taxon = "species",
                  min_age = "year",
                  max_age = "year",
                  value = "flagged")
  
  #Exclude problematic records
  dat_cl <- dat_cl[flags,]
  
  #Remove records with large uncertainty
  if ("coordinateUncertaintyInMeters" %in% colnames(dat_cl)) {
    dat_cl <- 
      dat_cl %>%
      filter(coordinateUncertaintyInMeters / 1000 <= 100 | is.na(coordinateUncertaintyInMeters))
  } else {
    warning("coordinateUncertaintyInMeters not found for ", unique(dat_cl$species))
  }
  
  #Remove unsuitable data data_sources, especially fossils 
  table(dat$basisOfRecord)
  dat_cl <- filter(dat_cl, basisOfRecord == "HUMAN_OBSERVATION" |
                     basisOfRecord == "OBSERVATION" |
                     basisOfRecord == "PRESERVED_SPECIMEN")
  
  # remove records before 2000
  dat_cl <- dat_cl %>% filter(year > 2000)
  
  return(dat_cl)
}

spp_search <- function(start, spp, stateProvince){
  occ <- occ_search(taxonKey = rgbif::name_backbone(spp)$usageKey, 
             hasCoordinate = T, 
             start = start, 
             stateProvince = stateProvince,
             limit = 1000)$data
  return(occ)
}
