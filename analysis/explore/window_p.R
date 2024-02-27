
#' Create moving window maps of allele frequencies
#'
#' @param x either a `vcfR` type object, a path to a .vcf file, or a allele frequency matrix
#' @param coords two-column matrix or data.frame representing x (longitude) and y (latitude) coordinates of samples
#' @param lyr RasterLayer to slide the window across
#' @param ... additional arguments to pass to `window_general`
#'
#' @return a RasterStack where each layer is a moving window allele frequency map
#' @export
#'
#' @examples
window_p2 <- function(x, coords, lyr, wdim = 3, fact = 0, ...){
  # convert to allele frequency matrix
  if (inherits(x, "vcfR") | is.character(x)) x <- wingen::vcf_to_dosage(x)/2

  # get moving window allele frequency for each allele
  if (!is.null(dim(x))) {
    r <-
      furrr::future_map(
        1:ncol(x),
        ~raster::raster(window_general(x = x[,.x], coords = coords, lyr = lyr, stat = mean, wdim = wdim, fact = fact, na.rm = TRUE,...)[[1]]),
        .options = furrr::furrr_options(
          seed = TRUE,
          packages = c("wingen", "terra", "raster")
        ),
        .progress = TRUE
      ) %>% 
      raster::stack() %>% 
      terra::rast()   
  } else {
    r <- wingen::window_general(x, coords = coords, lyr = lyr, stat = mean, wdim = wdim, fact = fact,  na.rm = TRUE, ...)
  }

  # give loci names to raster
  if (!is.null(colnames(x))) names(r) <- colnames(x)

  return(r)
}


window_p <- function(x, coords, lyr, wdim = 3, fact = 0, ...){
  # convert to allele frequency matrix
  if (inherits(x, "vcfR") | is.character(x)) x <- wingen::vcf_to_dosage(x)/2

  if (!is.null(dim(x))) {

    # Perform parallel computation using foreach
    r <- foreach(i = 1:ncol(x), .packages = c("wingen", "terra", "raster")) %dopar% {
      raster::raster(window_general(x = x[,i], coords = coords, lyr = lyr, stat = mean, wdim = wdim, fact = fact, na.rm = TRUE, ...)[[1]])
    } %>% 
    raster::stack() %>% 
    terra::rast()

  } else {
    r <- wingen::window_general(x, coords = coords, lyr = lyr, stat = mean, wdim = wdim, fact = fact,  na.rm = TRUE, ...)
  }

  # give loci names to raster
  if (!is.null(colnames(x))) names(r) <- colnames(x)

  return(r)
}


# function to run window_p and window_gd for each variable
run_windows <- function(snps_i, freqs, coords_proj, lyr){
  # create moving window map of allele freqs for each snp
  pstk <- window_p(freqs[,snps_i], coords = coords_proj, lyr = lyr, wdim = 11, fact = 0)

  # calculate moving window map of pi across all SNPs
  dpg <- window_general(freqs[,snps_i], coords = coords_proj, lyr = lyr, stat = "pi", wdim = 11, fact = 0, rarify = FALSE)

  return(list(pstk = pstk, dpg = dpg))
}

