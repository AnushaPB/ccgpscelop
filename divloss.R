
#' Quantify distinctness of genetic diversity across space using moving window maps of allele frequency
#'
#' @param x either a `vcfR` type object, a path to a .vcf file, or a dosage matrix
#' @param coords two-column matrix or data.frame representing x (longitude) and y (latitude) coordinates of samples
#' @param lyr RasterLayer to slide the window across
#' @param prop whether to return the proportional change (`TRUE`, default) or the raw change (`FALSE`)
#' @param binary whether to use binary presence/absence for alleles in each cell (`TRUE`) or use allele frequencies in each cell (`FALSE`, default). If `TRUE`, all cells in the moving window raster where the allele frequency is greater than 0 will be assigned a value of 1.
#' @param lower_cutoff allele frequency lower cutoff. Only alleles with a frequency above this value will be retained.
#' @param upper_cutoff allele frequency upper cutoff. Only alleles with a frequency below this value will be retained.
#' @param rare_allele whether to transform the genetic data such that the frequencies always reflect that of the rare allele at the locus (i.e., if p > 0.50, a 1 - p transformation would be applied)
#' @param ... additional arguments to pass to `window_general`
#'
#' @return
#' @export
#'
#' @examples
distinct_do_everything <- function(x, coords, lyr, lower_cutoff = NULL, upper_cutoff = NULL, prop = TRUE, binary = FALSE, rare_allele = FALSE, ...){
  # convert to dosage matrix
  if(inherits(x, "vcfR") | is.character(x)) x <- wingen::vcf_to_dosage(x)

  # generate stack of allele frequency rasters
  pstk <- window_p(x, coords, lyr, ...)

  # create distinct map for each layer
  # rare_only set to FALSE because transformation applied earlier
  dls <- distinct_map(pstk, lower_cutoff = lower_cutoff, upper_cutoff = upper_cutoff, prop = prop, binary = binary, rare_allele = rare_allele)

  return(c(dls, list(AlleleStack = pstk)))
}


#' Calculate distinct genetic diversity from allele frequency rasters
#'
#' @param x RasterStack of allele frequency maps produced by `window_p`
#' @param prop whether to return the proportional loss (TRUE, default) or the raw loss (FALSE)
#'
#' @return
#' @export
#'
#' @examples
distinct_map <- function(x, lower_cutoff = NULL, upper_cutoff = NULL, prop = TRUE, binary = FALSE, rare_allele = FALSE){
  # Transform allele frequency stack so the allele of interest (1) is always the minor/rarer allele
  if(rare_allele) x <- rare_pstk(x)

  # Filter with cutoffs and apply binary transformation
  if (!is.null(lower_cutoff) | !is.null(upper_cutoff) | binary) x <- transform_p(x, lower_cutoff, upper_cutoff, binary)

  # Create distinct map (weighted allele frequency map)
  dstk <- 
    map(terra::as.list(x), ~{
      p0 <- as.numeric(global(.x, fun = "sum", na.rm = TRUE))
      return(.x/p0)
    }, .progress = TRUE) %>%
    rast()

  # calculate the average loss across the entire stack
  davg <- mean(dstk, na.rm = TRUE)

  return(list(DistinctMean = davg,
              DistinctStack = dstk))
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

rare_pstk <- function(x){
  m <- unlist(global(x, fun = "mean", na.rm = TRUE))
  rare_only <- map(1:nlyr(x), ~if(m[.x] > 0.5) return(1 - x[[.x]]) else return(x[[.x]])) %>% rast()
  return(rare_only)
}

