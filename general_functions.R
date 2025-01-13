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

gglm <- function(x, y, df){
  ggplot(df, aes(x = .data[[x]], y = .data[[y]])) +
    geom_point(aes(col = .data[[x]])) +
    geom_smooth(method = "lm", col = "black") +
    labs(x = make_pretty_names(x), y = make_pretty_names(y)) +
    ggpubr::stat_cor(method = "pearson", label.y = min(df[[y]], na.rm = TRUE), label.x = min(df[[x]], na.rm = TRUE)) +
    theme_classic() +
    scale_color_viridis_c(option = "mako") +
    theme(legend.position = "none") 
}

ggpartial <- function(x, y, f, df){
  diff <- setdiff(f, x)
  response_f <- paste0(y, " ~ ", paste(diff, collapse = " + "))
  predictor_f <- paste0(x, " ~ ", paste(diff, collapse = " + "))
  
  response_resid <- residuals(lm(response_f, data = df))
  predictor_resid <- residuals(lm(predictor_f, data = df))

  df$response_resid <- response_resid
  df$predictor_resid <- predictor_resid

  ggplot(df, aes(x = predictor_resid, y = response_resid)) +
    geom_point(aes(col = .data[[x]])) +
    geom_smooth(method = "lm", col = "black") +
    labs(x = paste("Partial", make_pretty_names(x)), y = paste("Partial", y)) +
    ggpubr::stat_cor(method = "pearson", label.y = min(df$response_resid, na.rm = TRUE), label.x = min(df$predictor_resid, na.rm = TRUE)) +
    scale_color_viridis_c(option = "mako") +
    theme_classic() +
    theme(legend.position = "none")
}

make_pretty_names <- function(vars){
  map_chr(vars, \(x){
    if (x == "paleo_change") return("Paleoclimate temperature change (LIG)")
    if (x == "paleo_change_cur_lgm") return("Paleoclimate temperature change (LGM)")
    if (grepl("CHELSA_bio12", x)) return("Contemporary precipitation")
    if (grepl("CHELSA_bio1", x)) return("Contemporary temperature")
    if (grepl("csi", x)) return("Past climate stability")
    if (grepl("gHM", x)) return("Human modification")
    if (grepl("glacier", x)) return("Glacier (inside/outside)")
    if (grepl("Q", x)) return("Admixture")
    if (grepl("tmean_dif", x)) return("Contemporary temperature change")
    return(x)
  })
}
