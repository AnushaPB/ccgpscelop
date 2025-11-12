get_ca <- function() {
  # Load the U.S. state boundaries data
  states <- tigris::states(cb = TRUE)

  # Extract the boundary of California (CA)
  ca <- states[states$STUSPS == "CA", "STUSPS"]
  #ca <- st_read(here("data", "ca_state", "CA_State.shp"))
  #ca <- sf::st_transform(ca, sf::st_crs(4326))

  # Remove islands by taking the largest polygon
  ca_parts <- st_cast(ca, "POLYGON")
  ca_parts$area <- st_area(ca_parts)
  ca_mainland <- ca_parts[which.max(ca_parts$area), ]

  return(ca_mainland)
}

get_corrected_coords <- function() {
  cc <-
    readxl::read_excel(here("data", "CCGP_SCOC_coordfixes.xlsx")) %>%
    dplyr::select(SampleID = SequenceID, x = Longitude1, y = Latitude1)
}

get_coords <- function(sf = FALSE, all = FALSE) {
  # sample coords
  coords <- read_table(here("data", "ccgp_data", "58-Sceloporus.coords.txt"), col_names = FALSE)
  colnames(coords) <- c("SampleID", "x", "y")

  # CORRECTED COORDS:
  cc <- get_corrected_coords()

  # Add corrected coordinates
  coords <-
    coords %>%
    left_join(cc, by = "SampleID") %>%
    mutate(
      x = ifelse(is.na(x.y), x.x, x.y),
      y = ifelse(is.na(y.y), y.x, y.y)
    ) %>%
    dplyr::select(SampleID, x, y)

  # Check that coordinates were correctly replaced
  stopifnot(coords %>% filter(SampleID == "Sceocc_HBS142159") %>% pull(y) == cc %>% filter(SampleID == "Sceocc_HBS142159") %>% pull(y))

  # S. beckii samples
  beckii_samples <-
    c(
      "Scelocci_CCGPMC_MW01-3-14",
      "Scebec_7687",
      "Scebec_7727",
      "Scebec_7457",
      "Scebec_7499",
      "Scebec_7548",
      "Scebec_7553"
    )

  # Unknown provenance
  unknown_samples <- "Scelocci_CHI1382_DAW5-46-21"

  # Potentially swapped samples
  swapped_samples <- c("Scelocci_CAS213197", "Scelocci_CAS214858")

  # Filter for samples in the VCF
  fam <- read_table(here("data", "ccgp_data", "58-Sceloporus_complete_coords_annotated.fam"), col_names = FALSE)

  # Identify samples not in vcf
  ndropped <- length(setdiff(coords$SampleID, fam$X2))

  # Filter to samples in vcf
  coords <- coords %>% filter(SampleID %in% fam$X2)

  # If all = TRUE, return coordinates before removing beckii samples, unknown provenance sample, and potentially swapped samples
  if (all) {return(coords)}

  message(
    "Removing: ", length(beckii_samples), " S. beckii samples, ",
    length(unknown_samples), " unknown provenance sample, and ",
    length(swapped_samples), " potentially swapped samples",
    "\nNumber of samples dropped from VCF during QC: ", ndropped
  )

  # Remove beckii, unknown provenance sample, and potentially swapped samples
  coords <-
    coords %>%
    filter(!SampleID %in% c(beckii_samples, unknown_samples, swapped_samples))

  if (sf) {
    coords <- st_as_sf(coords, coords = c("x", "y"), crs = 4326)
  }

  message(nrow(coords), " samples with coordinates")

  return(coords)
}

get_range <- function(){
  range_map <-
    st_read(here("data", "rWFLIx_CONUS_HabMap_2001v1", "rWFLIx_CONUS_Range_2001v1.shp")) %>%
    st_transform(3310)# %>%
    # st_intersection(get_ca() %>% st_transform(3310))
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

get_pops <- function(){
  path <- here("analysis", "admixture", "outputs", "Q8.csv")
  message("Reading in ", path)
  pop_df <- read_csv(path) %>% distinct(cluster, SampleID) %>% mutate(cluster = factor(cluster))
}

get_pop_cols <- function(){
  cs1 <- viridis::mako(4, begin = 0.3, end = 0.8)
  names(cs1) <- c("7", "5", "4", "1")
  cs2 <- viridis::rocket(9-4, begin = 0.3, end = 0.9)
  names(cs2) <- c("9", "8", "6", "2", "3")
  cs <- c(cs1, cs2)
  return(cs)
}

gglm <- function(x, y, df, col = NULL){
  if (is.null(col)) col <- x
  ggplot(df, aes(x = .data[[x]], y = .data[[y]])) +
    geom_point(aes(col = .data[[col]])) +
    geom_smooth(method = "lm", col = "black") +
    labs(x = make_pretty_names(x), y = make_pretty_names(y)) +
    ggpubr::stat_cor(method = "pearson", label.y = min(df[[y]], na.rm = TRUE), label.x = min(df[[x]], na.rm = TRUE)) +
    theme_classic() +
    scale_color_viridis_c(option = "mako") +
    theme(legend.position = "none")
}

ggpartial <- function(x, y, f, df, col = NULL, alpha = 1, cex = 1){
  diff <- setdiff(f, x)
  response_f <- paste0(y, " ~ ", paste(diff, collapse = " + "))
  predictor_f <- paste0(x, " ~ ", paste(diff, collapse = " + "))

  response_resid <- residuals(lm(response_f, data = df))
  predictor_resid <- residuals(lm(predictor_f, data = df))

  df$response_resid <- response_resid
  df$predictor_resid <- predictor_resid

  if (!is.null(col)) {
    df$color <- df[[col]]
  } else {
    df$color <- 0
  }

  pretty_name_x <- make_pretty_names(x)
  pretty_name_x_lower <- paste0(tolower(substr(pretty_name_x, 1, 1)), substr(pretty_name_x, 2, nchar(pretty_name_x)))
  if (pretty_name_x == "NDVI") pretty_name_x_lower <- "NDVI"

  pretty_name_y <- make_pretty_names(y)
  pretty_name_y_lower <- paste0(tolower(substr(pretty_name_y, 1, 1)), substr(pretty_name_y, 2, nchar(pretty_name_y)))


  plt <-
    ggplot(df, aes(x = predictor_resid, y = response_resid)) +
    #geom_point(aes(fill = color), pch = 21, alpha = alpha, cex = cex, col = NA) +
    geom_point(aes(col = color), alpha = alpha, cex = cex) +
    geom_smooth(method = "lm", col = "black") +
    labs(x = paste("Partial", pretty_name_x_lower), y = paste("Partial", pretty_name_y_lower)) +
    ggpubr::stat_cor(method = "pearson", label.y = min(df$response_resid, na.rm = TRUE), label.x = min(df$predictor_resid, na.rm = TRUE)) +
    theme_classic() +
    theme(legend.position = "none")

  if (!is.null(col)) {
    plt <- plt + scale_color_viridis_c(option = "mako")
  } else {
    plt <- plt + scale_color_gradient(low = "darkgray", high = "darkgray")
  }

  return(plt)
}


ggpartialsem <- function(x, y, f, df, listw, col = TRUE, alpha = 1, cex = 1){
  diff <- setdiff(f, x)
  response_f <- paste0(y, " ~ ", paste(diff, collapse = " + "))
  predictor_f <- paste0(x, " ~ ", paste(diff, collapse = " + "))

  df$response_resid <- residuals(errorsarlm(response_f, listw = listw, data = df, zero.policy = TRUE))
  df$predictor_resid <- residuals(errorsarlm(predictor_f, listw = listw, data = df, zero.policy = TRUE))

  ggplot(df, aes(x = predictor_resid, y = response_resid)) +
    geom_point(aes(col = .data[[x]]), alpha = alpha, cex = cex) +
    geom_smooth(method = "lm", col = "black") +
    labs(x = paste("Partial", make_pretty_names(x)), y = paste("Partial", make_pretty_names(y))) +
    ggpubr::stat_cor(method = "pearson", label.y = min(df$response_resid, na.rm = TRUE), label.x = min(df$predictor_resid, na.rm = TRUE)) +
    scale_color_viridis_c(option = "mako") +
    theme_classic() +
    theme(legend.position = "none")
}

ggmap <- function(x){
  x_df <- as.data.frame(x, xy = TRUE)
  ca <- get_ca()
  ggplot(x_df) +
    geom_sf(data = ca, fill = "white", color = "black") +
    geom_raster(aes(x = x, y = y, fill = .data[[names(x)]])) +
    labs(fill = make_pretty_names(names(x))) +
    theme_void() +
    scale_fill_viridis_c(option = "mako")
}

make_pretty_names <- function(vars){
  map_chr(vars, \(x){
    if (x == "lgm") return("LGM")
    if (x == "lig") return("LIG")
    if (x == "cur_lig") return("Paleoclimate change\n(CUR - LIG)")
    if (x == "cur_lgm") return("Paleoclimate change\n(CUR - LGM)")
    if (x == "lgm_lig") return("Paleoclimate change\n(LGM - LIG)")
    if (x == "eh_lgm") return("Paleoclimate change\n(EH - LGM)")
    if (x == "csi_custom") return("Temperature stability\n(MIS19 to CUR)")
    if (grepl("bio1", x)) return("Contemporary temperature")
    if (grepl("csi_past", x)) return("Paleoclimate stability")
    if (grepl("csi_future", x)) return("Future climate stability")
    if (grepl("gHM", x)) return("Human modification")
    if (grepl("glacier", x)) return("Glacier (inside/outside)")
    if (grepl("Q", x)) return("Admixture")
    if (grepl("tmean_dif", x)) return("Contemporary temperature change")
    if (grepl("fire", x)) return("Wildland fire frequency")
    if (grepl("resid", x)) return("Residuals")
    if (grepl("genes_ho", x) | grepl("Ho_gea", x)| grepl("gea_Ho", x)) return("Adaptive heterozygosity")
    if (grepl("Ho_syn", x)) return("Synonymous heterozygosity")
    if (grepl("Ho_nonsyn", x)) return("Non-synonymous heterozygosity")
    if (grepl("Ho", x)) return("Genome-wide heterozygosity")
    return(x)
  })
}

scaffold_theme <- function() {
  list(
    scale_x_continuous(labels = scales::comma, expand = c(0, 0)),
    scale_y_continuous(labels = scales::comma),
    theme_classic(),
    facet_wrap(~scaffold, scales = "fixed", ncol = 1),
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      strip.background = element_blank(),
      strip.text.x = element_text(hjust = 0, size = 10),
      legend.position = "bottom",
    )
  )
}


scaffold_theme_y <- function() {
  list(
    scale_x_continuous(labels = scales::comma, expand = c(0, 0)),
    scale_y_continuous(labels = scales::comma),
    theme_classic(),
    facet_wrap(~scaffold, scales = "fixed", ncol = 1),
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.line.x = element_blank(),
      strip.background = element_blank(),
      strip.text.x = element_text(hjust = 0, size = 10),
      legend.position = "bottom",
    )
  )
}

#' Import environmental layers of choice for RDA, adaptive index, or genomic offset calculations
#'
#' @param type options are "pca" or "bio1ndvi" for BIO1 + NDVI
#' @param future whether to also import future env layers (defaults to FALSE)
#'
#' @return
#' @export
get_envlayers <- function(type = "bio1ndvi", future = FALSE) {
  if (type == "pca") {
    env_pres <- terra::rast(here("data", "env", "california_chelsa_bioclim_1981-2010_V.2.1_pca.tif"))
      # terra::project("epsg:3310")
    names(env_pres) <- paste("env_", names(env_pres), sep = "")
  }
  if (type == "bio1ndvi") {
    bioclim <- rast(here("data", "env", "california_chelsa_bioclim_1981-2010_V.2.1.tif"))
    #  %>% project("epsg:3310")
    bio1 <- bioclim[["CHELSA_bio1_1981-2010_V.2.1"]]
    ndvi <- terra::rast(here("data", "env", "california_ndvi_mean_2000_2020.tif"))
    ndvi <- terra::resample(ndvi, bio1, method = "bilinear")
    env_pres <- c(bio1, ndvi)
    names(env_pres) <- c("BIO1", "NDVI")
  }

  if (future) {
    if (type == "pca") {
      env_fut_1 <- terra::rast(paste0(here("data", "env", "future"), "/CHELSA_2071-2100_", cap_model, "_", ssp[1], "_V.2.1_pca.tif"))
      env_fut_2 <- terra::rast(paste0(here("data", "env", "future"), "/CHELSA_2071-2100_", cap_model, "_", ssp[2], "_V.2.1_pca.tif"))
      env_fut <- c(env_fut_1, env_fut_2)
    }

    if (type == "bio1ndvi") {
      env_fut <- terra::rast(here("data", "env", "future", "env_fut_2071-2100_GFDL-ESM4_ssp126_ssp585.tif"))
      names(env_fut) <- c("CHELSA_bio1_2071-2100_gfdl-esm4_ssp126_V.2.1", "CHELSA_bio1_2071-2100_gfdl-esm4_ssp585_V.2.1", "NDVI")
    }
  }
  if (!future) env_fut <- NULL
  return(list(env_pres = env_pres, env_fut = env_fut))
}
