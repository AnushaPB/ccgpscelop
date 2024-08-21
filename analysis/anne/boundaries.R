# Structure boundaries ----------------------------------------------------
# -------------------------------------------------------------------------

#' Make admixture zone boundaries for each K value
#'
#' @param Kvalue K to calculate
#' @param gg_df tidy dataframe of admixture results
#' @param threshold threshold for admixture zone; if set to NULL, all Q-values will be mapped
#' @param grid raster for resampling
#'
#' @return
#' @export
make_boundary <- function(Kvalue, krig_admix, threshold = c(0.2, 0.6), grid) {
  rs <-  krig_admix[[1]]
  rs[] <- NA

  # Set up tidy df
  admix <- krig_admix %>%
    terra::as.data.frame(x, xy = TRUE, na.rm = FALSE) %>%
    tidyr::as_tibble() %>%
    tidyr::pivot_longer(names_to = "K", values_to = "Q", -c(x, y)) %>%
    dplyr::mutate(K = as.factor(gsub("K", "", K))) %>%
    dplyr::group_by(x, y) %>%
    drop_na(Q) %>% # remove NAs
    dplyr::group_by(x, y) %>% # should already be grouped but just to be safe
    dplyr::filter(K == Kvalue) %>%
    dplyr::filter(Q > threshold[1] & Q < threshold[2]) # to just get admixture Q-values
  # dplyr::mutate(newQ = case_when(Q > threshold[1] & Q < threshold[2] ~ Q,
  #                                Q > threshold[2] ~ Kvalue))

  # layer 1 is just K value; layer 2 are Q-values; layer 3 is newQ
  admixrl <- terra::rast(admix, type = "xyz")

  # Ensure extents will be the same for mosaic-ing later
  if (!inherits(grid, "SpatRaster")) grid <- terra::rast(grid)
  admixrl <- terra::resample(admixrl, grid)

  return(admixrl[[1]])
}

#' Convert genoscape raster stack to polygons
#' From https://github.com/mgdesaix/mignette/blob/main/R/scape_to_shape.R
#'
#' @param x Genoscape RasterStack
#' @param prob_threshold Probability value to include raster cell in polygon
#' @param d Distance threshold for smoothr::drop_crumbs()
#' @param f Distance threshold for smoothr::fill_holes()
#' @param s Smoothness
#'
#' @return A polygon sf object for all genoscape clusters
#' @export
#'
scape_to_shape <- function(x, prob_threshold = 0.5, d = 1000, f = 1000, s = 3){
  m <- c(-Inf, prob_threshold, NA)
  rclmat <- matrix(m, ncol = 3, byrow = TRUE)
  sf::sf_use_s2(FALSE)
  genoscape_classified <- terra::classify(x, rclmat, right = FALSE) %>%
    terra::as.list()

  rtp <- function(y){
    poly.tmp <- terra::as.polygons(y)
    poly.tmp.cluster <- names(poly.tmp)
    names(poly.tmp) <- "Cluster"
    poly.tmp[[1]] <- poly.tmp.cluster
    poly.tmp <- terra::project(poly.tmp, y="+proj=longlat +datum=WGS84")
    poly.tmp <- sf::st_as_sf(poly.tmp)

    poly.tmp.dc <- smoothr::drop_crumbs(poly.tmp,
                                        threshold = units::set_units(d, km^2))
    poly.tmp.fh <- smoothr::fill_holes(poly.tmp.dc,
                                       threshold = units::set_units(f,km^2))
    poly.tmp.smooth <- smoothr::smooth(poly.tmp.fh, method = "ksmooth", smoothness = s)
    poly.tmp.smooth <- sf::st_as_sf(poly.tmp.smooth, 'Spatial')
    return(poly.tmp.smooth)
  }
  polygon.list <- lapply(genoscape_classified, rtp)
  polygon.sf <- do.call("rbind", polygon.list)
  return(polygon.sf)
}
