
get_het <- function(){
  # get # callable sites for denominator
  callable <- read.csv(here("data", "ccgp_data", "callable_counts.csv"))
  blra_callable <- callable %>% filter(species == "Laterallus jamaicensis") %>% pull(callable_sites_post_filter)
  vira_callable <- callable %>% filter(species == "Rallus limicola") %>% pull(callable_sites_post_filter)

  het_blra <- format_het(here("analysis", "genetic_diversity", "outputs", "blra.het"), blra_callable)
  het_blra$SampleID <- recode_blra(het_blra$IID)
  het_blra$GenusSpecies <- "Laterallus jamaicensis"

  het_vira <- format_het(here("analysis", "genetic_diversity", "outputs", "vira.het"), vira_callable)
  het_vira$SampleID <- recode_vira(het_vira$IID)
  het_vira$GenusSpecies <- "Rallus limicola"

  return(bind_rows(het_blra, het_vira))
}

format_het <- function(path, callable_sites){
  # Load the data in R
  het_data <- read.table(path, header = TRUE)

  # Calculate the average heterozygosity per individual
  het_data$Ho <- (het_data$N.NM. - het_data$O.HOM.)/ callable_sites

  return(het_data)
}

gd_mod <- function(df, y = "Ho"){
  # make site df
  sub_sites <- 
    df %>%
    distinct(x, y, .keep_all = TRUE) %>%
    mutate(site = 1:nrow(.)) %>%
    st_as_sf(coords = c("x", "y"))
  
  # calculate gravity model metrics
  gm_sf <- gm_metrics(sub_sites)
  gm_df <- data.frame(as_Spatial(gm_sf))
    
  # dfine variables 
  cat_cols <- c("wetland_type", "wetland_subtype", "vegetation_modifier", "tidal_modifier")
  cont_cols <- c("degree", "closeness", "Shape_Area", "Annual.Precipitation", "betweenness")

  # transform data
  mod_df <- 
    st_join(df, sub_sites) %>%
    as_Spatial() %>%
    data.frame() %>%
    left_join(gm_df) %>%
    mutate_at(all_of(cont_cols), scale)

  # check if there are at least two groups with more than one sample in each group OR if there are more than 2 groups
  good_vars <- which(map_lgl(mod_df[, cat_cols], ~sum((table(.x) > 1)) > 2 | (length(unique(.x)) > 2)))
  cat_cols <- cat_cols[good_vars]

  # make model formula 
  f <- formula(paste(y, "~", paste(c(cont_cols, cat_cols), collapse = "+")))
  mod <- data.frame(summary(lm(f, data = mod_df))$coefficients)
  mod$Variable <- row.names(mod)

  out <- data.frame(mod, gd = y)
  return(out)
}


visualize_vars <- function(df, vars, model = NULL, GenusSpecies = NULL, group = NULL){

  if (is.null(GenusSpecies) & is.null(group)){
    GenusSpecies <- ifelse(grepl("BLRA", model), "Laterallus jamaicensis", "Rallus limicola")
    group <- ifelse(grepl("Bay Area", model), "Bay Area", "Sierra Foothills")
  }

  # make site df
  sub_sites <- 
    df %>%
    filter(GenusSpecies == {{GenusSpecies}}, group == {{group}}) %>%
    distinct(x, y, .keep_all = TRUE) %>%
    mutate(site = 1:nrow(.)) %>%
    st_as_sf(coords = c("x", "y"))
  
  # calculate gravity model metrics
  gm_sf <- gm_metrics(sub_sites)
  # make network graph
  gg <- graph2nb(gabrielneigh(st_coordinates(sub_sites)), sym = TRUE)
  gg <- nb2lines(gg, coords = sf::st_coordinates(sub_sites), proj4string = st_crs(sub_sites), as_sf=TRUE)

  # get variables
  gm_gg <- 
    gm_sf %>%
     dplyr::select(all_of(vars)) %>%
     mutate_at(vars, scale) %>%
     pivot_longer(all_of(vars)) %>%
     drop_na(value)

  # make background
  ca <- get_ca()
  sa <- st_crop(ca, st_bbox(st_buffer(gm_gg, 1000)))

  plt <- 
   ggplot(gm_gg) +
   geom_sf(data = sa) +
   geom_sf(data = gg, col = "darkgray") +
   geom_sf(aes(fill = value), col = "black", pch = 21, cex = 3) +
   facet_wrap(~name, strip.position = "bottom") +
   theme_void() +
   scale_fill_viridis_c(option = "plasma") +
   theme(
    legend.position = "none", 
    strip.text=element_text(size = 12, vjust = -0.1),
    plot.title = element_text(hjust = 0.5, face = "italic")) +
   ggtitle(paste0({{GenusSpecies}}, " (", {{group}}, ")"))

  return(plt)
}


gd_mod2 <- function(df){
  # make site df
  sub_sites <- 
    df %>%
    distinct(x, y, .keep_all = TRUE) %>%
    mutate(site = 1:nrow(.)) %>%
    st_as_sf(coords = c("x", "y"))
  
  # calculate gravity model metrics
  gm_sf <- gm_metrics(sub_sites)
  gm_df <- data.frame(as_Spatial(gm_sf))
    
  # dfine variables 
  cat_cols <- c("wetland_type", "wetland_subtype", "vegetation_modifier", "tidal_modifier")
  cont_cols <- c("degree", "closeness", "Shape_Area", "Annual.Precipitation", "betweenness")

  # make model df
  mod_df <- 
    st_join(df, sub_sites) %>%
    as_Spatial() %>%
    data.frame() %>%
    left_join(gm_df) %>%
    mutate_at(all_of(cont_cols), scale)
  
  # check if there are at least two groups with more than one sample in each group OR if there are more than 2 groups
  good_vars <- which(map_lgl(mod_df[, cat_cols], ~sum((table(.x) > 1)) > 2 | (length(unique(.x)) > 2)))
  cat_cols <- cat_cols[good_vars]

  interactions <- paste0("group*", c(cont_cols, cat_cols))
  f <- formula(paste("Ho ~", paste(interactions, collapse = "+")))
  Ho_mod <- data.frame(summary(lm(f, data = mod_df))$coefficients)
  Ho_mod$Variable <- row.names(Ho_mod)

  out <- data.frame(Ho_mod, gd = "Ho")
  return(out)
}


get_roh <- function(){
  sub_coords <- get_coords(filter = TRUE)

  blra_ba_roh <- 
    read.table(here("analysis", "genetic_diversity", "outputs", "blra_ba.froh")) %>%
    mutate(SampleID = recode_blra(V1), froh = V2) %>%
    right_join(filter(sub_coords$blra, group == "Bay Area")) %>%
    # missing froh values indicates no roh greater than the minimum size were found
    mutate(froh0 = case_when(is.na(froh) ~ 0, .default = froh)) %>%
    mutate(spp = "Laterallus jamaicensis")

  blra_sf_roh <- 
    read.table(here("analysis", "genetic_diversity", "outputs", "blra_sf.froh")) %>%
    mutate(SampleID = recode_blra(V1), froh = V2) %>%
    right_join(filter(sub_coords$blra, group == "Sierra Foothills")) %>%
    mutate(froh0 = case_when(is.na(froh) ~ 0, .default = froh)) %>%
    mutate(spp = "Laterallus jamaicensis")

  vira_roh <- 
    read.table(here("analysis", "genetic_diversity", "outputs", "vira.froh")) %>%
    mutate(SampleID = recode_vira(V1), froh = V2) %>%
    right_join(filter(sub_coords$vira, group == "Sierra Foothills")) %>%
    mutate(froh0 = case_when(is.na(froh) ~ 0, .default = froh)) %>%
    mutate(spp = "Rallus limicola")

  roh_df <- 
    bind_rows(vira_roh, blra_ba_roh, blra_sf_roh) %>% 
    mutate(group = factor(group, levels = c("Sierra Foothills", "Bay Area"))) %>%
    dplyr::select(-V1, -V2) %>%
    st_as_sf()
    
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
    group_by(x, y, GenusSpecies, group) %>%
    summarize_at({{stat}}, mean, na.rm = TRUE) %>%
    drop_na(x, y) %>%
    st_as_sf(coords = c("x", "y"), crs = st_crs(coords)) %>%
    st_transform(4326)

}

# pop gen stats ---------------

het.dosage <- function(x, L = NULL){
  y <- x
  y[x == 0] <- 0
  y[x == 2] <- 0
  y[x == 1] <- 1

if (is.null(dim(y))) {
    if (!is.null(L)) y <- c(y, rep(0, L - length(x)))
    gd <- mean(y, na.rm = TRUE)
  } else {
    het_by_locus <- colMeans(y, na.rm = TRUE)
    if (!is.null(L)) het_by_locus <- c(het_by_locus, rep(0, L - ncol(x)))
    gd <- mean(het_by_locus, na.rm = TRUE)
  }

  return(gd)
}

#' Calculate mean heterozygosity
#'
#' @param hetmat matrix of heterozygosity (0/FALSE = homozygote, 1/TRUE = heterozygote)
#'
#' @return heterozygosity averaged across all individuals then all loci
#'
#' @noRd
calc_mean_het <- function(x) {

  if (inherits(x, "vcfR")) x <- vcf_to_het(x)

  if (is.null(dim(x))) {
    gd <- mean(x, na.rm = TRUE)
  } else {
    het_by_locus <- colMeans(x, na.rm = TRUE)
    gd <- mean(het_by_locus, na.rm = TRUE)
  }
  names(gd) <- "Ho"
  return(gd)
}

#' Convert vcf to heterozygosity matrix
#'
#' @param x can either be an object of class 'vcfR' or a path to a .vcf file
#'
#' @return heterozygosity matrix
#'
#' @noRd
vcf_to_het <- function(x) {
  het <- vcfR::is.het(vcfR::extract.gt(x), na_is_false = FALSE)

  # IMPORTANT: transform matrix so that rows are individuals and cols are loci
  het <- t(het)

  # if gen is a vector of only one site, turn into matrix with one column
  if (nrow(x@gt) == 1) {
    het <- matrix(het, ncol = 1)
  }

  return(het)
}






