
#' Quantify potential diversity loss across space using moving window maps of allele frequency
#'
#' @param x either a `vcfR` type object, a path to a .vcf file, or a dosage matrix
#' @param coords two-column matrix or data.frame representing x (longitude) and y (latitude) coordinates of samples
#' @param lyr RasterLayer to slide the window across
#' @param prop whether to return the proportional change (`TRUE`, default) or the raw change (`FALSE`)
#' @param binary whether to use binary presence/absence for alleles in each cell (`TRUE`) or use allele frequencies in each cell (`FALSE`, default). If `TRUE`, all cells in the moving window raster where the allele frequency is greater than 0 will be assigned a value of 1.
#' @param lower_cutoff allele frequency lower cutoff. Only alleles with a frequency above this value will be retained.
#' @param upper_cutoff allele frequency upper cutoff. Only alleles with a frequency below this value will be retained.
#' @param ... additional arguments to pass to `window_general`
#'
#' @return
#' @export
#'
#' @examples
divloss <- function(x, coords, lyr, lower_cutoff = NULL, upper_cutoff = NULL, prop = TRUE, binary = FALSE, ...){

  # convert to dosage matrix
  if(inherits(x, "vcfR") | is.character(x)) x <- wingen::vcf_to_dosage(x)

  # make sure that the rarest allele is coded as 1
  x <- rare_dos(x)

  # generate stack of allele frequency rasters
  pstk <- window_p(x, coords, lyr, ...)

  # calculate diversity loss index for each layer
  dlstk <- divloss_p(pstk, lower_cutoff = lower_cutoff, upper_cutoff = upper_cutoff, prop = prop, binary = binary)

  # calculate the average loss across the entire stack
  dlavg <- mean(dlstk, na.rm = TRUE)

  return(list(AvgLoss = dlavg,
              LossRasterStack = dlstk,
              AlleleRasterStack = pstk))
}


#' Create diversity loss rasters from allele frequency rasters
#'
#' @param x RasterStack of allele frequency maps produced by `window_p`
#' @param prop whether to return the proportional loss (TRUE, default) or the raw loss (FALSE)
#'
#' @return
#' @export
#'
#' @examples
divloss_p <- function(x, lower_cutoff = NULL, upper_cutoff = NULL, prop = TRUE, binary = FALSE){

  # filter with cutoffs and apply binary transformation
  if(!is.null(lower_cutoff) | !is.null(upper_cutoff) | binary) x <- transform_p(x, lower_cutoff, upper_cutoff, binary)

  # wrap up for parallel task
  x <- terra::as.list(x)
  wx <- purrr::map(x, terra::wrap)

  # run distinctnesss calculations
  wr <- furrr::future_map(wx, ~divloss_p_helper(x = .x, prop = prop))

  # unwrap SpatRaster
  r <- purrr::map(wr, terra::unwrap)
  r <- terra::rast(r)

  return(r)
}


#' Transform allele frequency rasters
#'
#' @param x RasterStack of allele frequency maps produced by `window_p`
#' @param lower_cutoff allele frequency lower cutoff. Only alleles with a frequency above this value will be retained.
#' @param upper_cutoff allele frequency upper cutoff. Only alleles with a frequency below this value will be retained.
#' @param binary whether to use binary presence/absence for alleles in each cell (`TRUE`) or use allele frequencies in each cell (`FALSE`, default). If `TRUE`, all cells in the moving window raster where the allele frequency is greater than 0 will be assigned a value of 1.
#'
#' @return
#' @export
#'
#' @examples
transform_p <- function(x, lower_cutoff = NULL, upper_cutoff = NULL, binary = FALSE){
  # convert layers to binary
  if (binary) x[x > 0] <- 1

  # get mean (global p) from each layer
  p <- purrr::map(terra::as.list(x), ~global(.x, "mean", na.rm = TRUE)[1,1])

  # remove alleles based on cutoffs
  if(!is.null(upper_cutoff)) x <- x[[which(p < upper_cutoff)]]
  if(!is.null(lower_cutoff)) x <- x[[which(p > lower_cutoff)]]

  return(x)
}



#' Helper function for divloss_p
#' @param wx packed SpatRaster
#' @inheritParams divloss_p
#' @noRd
divloss_p_helper <- function(wx, lower_cutoff = NULL, upper_cutoff = NULL, prop = TRUE){
  # convert back to SpatRaster 
  x <- unwrap(wx)

  # Refill raster with values produced by calc_loss
  x[] <- purrr::map_dbl(1:terra::ncell(x), ~calc_loss(i = .x, x = x, prop = prop))

  # wrap up again
  x <- wrap(x)

  return(x)
}

#' Helper function to calculate loss statistic
#' @param i raster cell index
#' @param x RasterLayer of allele frequencies
#' @inheritParams divloss_p
#' @noRd
calc_loss <- function(i, x, prop = TRUE){

  # if the cell has an NA value, return NA
  if(is.na(x[i])) return(NA)

  # calculate the initial global allele frequency
  begin <- terra::global(x, "mean", na.rm = TRUE)

  # assign the cell a value of 0 (e.g. "remove" that cell)
  # note: you don't assign this cell an NA value because you
  # want the total number of cells used to calculate the mean
  # to remain the same (otherwise you could end up with positive change)
  x[i] <- 0

  # calculate the resulting global allele frequency
  end <- terra::global(x, "mean", na.rm = TRUE)

  # calculate the change
  change <- end - begin

  # if prop, then divide the change by the initial frequency
  if (prop) change <- change/begin

  # convert to dbl
  change <- change[1,1]

  return(change)
}

rare_dos <- function(dos) apply(dos, 2, rare_dos_helper)

rare_dos_helper <- function(x){
  nalt <- sum(x == 2)*2 + sum(x == 1)
  nref <- sum(x == 0)*2 + sum(x == 1)
  if (nalt > nref) return(2 - x) else return(x)
}


