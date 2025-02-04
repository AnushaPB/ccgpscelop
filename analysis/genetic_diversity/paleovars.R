get_paleovars <- function(paleovars = c("mis19", "lig", "lgm", "hs1", "ba", "yds", "eh", "mh", "lh", "cur"), cache = FALSE, process = TRUE) {
  ca <- get_ca()
  # download data
  walk(paleovars, ~ {
    rpaleoclim::paleoclim(.x, "2_5m", region = ext(ca), cache_path = here("data", "env", "paleoclim"))
  })

  # Define the file names and paths
  file_names <- c("MIS19_v1_r2_5m.zip", "LIG_v1_2_5m.zip", "chelsa_LGM_v1_2B_r2_5m.zip", "HS1_v1_2_5m.zip", "BA_v1_2_5m.zip", "YDS_v1_2_5m.zip", "EH_v1_2_5m.zip", "MH_v1_2_5m.zip", "LH_v1_2_5m.zip", "CHELSA_cur_V1_2B_r2_5m.zip")
  names(file_names) <- c("mis19", "lig", "lgm", "hs1", "ba", "yds", "eh", "mh", "lh", "cur")
  file_paths <- here("data", "env", "paleoclim", file_names[paleovars])

  # Load and process the data using purrr
  if (process) {
    env_list <- purrr::map(file_paths, ~ {
      data <- rpaleoclim::load_paleoclim(.x)
      data <- terra::crop(data, ca)
      data <- terra::mask(data, ca)
      data
  })
  } else {
    env_list <- purrr::map(file_paths, ~ {
      data <- rpaleoclim::load_paleoclim(.x)
      data
    })
  }
 
  names(env_list) <- paleovars

  # Cache data
  path =  here("data", "env", "paleoclim", "paleoclim.tif")
  if (!file.exists(path) && cache){
    # Add unique names
    # Note: doesn't work to just turn list into stack, the names will be the list names plus the index
    env_list_names <- unlist(imap(env_list, ~paste0(.y, "_", names(.x))))
    env_stack <- rast(env_list)
    names(env_stack) <- env_list_names
    writeRaster(env_stack, path, overwrite = TRUE)
  }

  return(env_list)
}

get_paleokey <- function(){
  # Create a key with layer information
  paleokey <- data.frame(
    id = c("cur", "lh", "mh", "eh", "yds", "ba", "hs1", "lgm", "lig", "mis19", "mpwp", "m2"),
    period = c("Current", "Late Holocene: Meghalayan", "Mid Holocene: Northgrippian", "Early Holocene: Greenlandian", "Pleistocene: Younger Dryas Stadial", "Pleistocene: Bølling-Allerød", "Pleistocene: Heinrich Stadial 1", "Pleistocene: Last Glacial Maximum", "Pleistocene: Last Interglacial", "Pleistocene: MIS19", "Pliocene: Mid-Pliocene warm period", "Pliocene: M2"),
    time = c("1979 – 2013", "4.2-0.3 ka", "8.326-4.2 ka", "11.7-8.326 ka", "12.9-11.7 ka", "14.7-12.9 ka", "17.0-14.7 ka", "ca. 21 ka", "ca. 130 ka", "ca. 787 ka", "3.205 Ma", "ca. 3.3 Ma"),
    source = c("CHELSA", "Fordham et al. 2017", "Fordham et al. 2017", "Fordham et al. 2017", "Fordham et al. 2017", "Fordham et al. 2017", "Fordham et al. 2017", "CHELSA", "Otto-Bliesner et al. 2006", "Brown et al. 2018", "Hill 2015", "Dolan et al. 2015")
  )
}
