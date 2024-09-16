
# Admixture functions -----------------------------------------------------
# -------------------------------------------------------------------------

#' Get CV from admixture run
#'
#' @param output_dir directory with admixture files
#' @param prefix file prefix
#'
#' @return
#' @export
get_cv <- function(output_dir, prefix = "58-Sceloporus_maf05_minDP5_maxDP50_rmsamp60_mm80_rmsamp20_r60"){
  # create a vector with file names
  file_names <- paste0(output_dir, paste0("/", prefix, ".", 1:10, ".out"))

  # get cv errors
  safe_read <- possibly(readLines)
  cv_errors <- purrr::map_dbl(file_names, ~ {
    file_content <- safe_read(.x)
    if (is.null(file_content)) return(NA)
    cv_line <- grep("CV error", file_content, value = TRUE)
    if (length(cv_line) > 0) {
      cv_value <- as.numeric(sub(".*CV error \\(K=[0-9.]+\\): ([0-9.]+).*", "\\1", cv_line))
      cv_value
    } else {
      NA
    }
  })

  df <- data.frame(K = factor(1:length(cv_errors), levels = 1:length(cv_errors)), cv_error = cv_errors)

  return(df)
}

#' Get Q-values for specific K
#'
#' @param K K-value
#' @param output_dir directory with admixture files
#' @param prefix file prefix
#' @param qmat_only if TRUE, only return Q-values; if FALSE, also return individual max assignment into cluster and sample ID
#'
#' @return
#' @export
get_Q <- function(K, output_dir, prefix = "58-Sceloporus_maf05_minDP5_maxDP50_rmsamp60_mm80_rmsamp20_r60", qmat_only = FALSE){
  # use fam to get sampleID order/names
  fam <- data.frame(read_table(paste0(output_dir, paste0("/", prefix, ".fam")), col_names = FALSE))

  # read in Q vals
  qmat <- read_table(paste0(output_dir, paste0("/", prefix, ".", K, ".Q")), col_names = FALSE)

  if (qmat_only) return(qmat)

  # add clusters
  qmat$cluster <- factor(apply(qmat, 1, which.max))

  # add SampleID
  qmat$SampleID <- as.vector(fam[ ,2])

  return(qmat)
}

#' Build structure-style bar plot for admixture results
#'
#' @param qmat df of Q values
#' @param sort_by_Q whether to sort individuals by max Q value
#' @param legend whether to plot legend
#'
#' @return
#' @export
structure_plot <- function(qmat, sort_by_Q = TRUE, legend = TRUE) {
  # Get K
  K <- ncol(qmat)

  dat <- as.data.frame(qmat)
  dat <- dat %>%
    tibble::rownames_to_column(var = "order")

  if (sort_by_Q) {
    gr <- apply(qmat, MARGIN = 1, which.max)
    gm <- max(gr)
    gr.o <- order(sapply(1:gm, FUN = function(g) mean(qmat[, g])))
    gr <- sapply(gr, FUN = function(i) gr.o[i])
    or <- order(gr)

    dat <- dat %>%
      dplyr::arrange(factor(order, levels = or))
    dat$order <- factor(dat$order, levels = dat$order)
  }

  # Make into tidy df
  gg_df <-
    dat %>%
    tidyr::pivot_longer(names_to = "cluster", values_to = "Q_value", -c(order)) %>%
    dplyr::mutate(cluster = gsub("X", "", cluster))

  # Build plot using helper function
  plt <- ggbarplot_helper(gg_df) + scale_fill_manual(values = viridis::turbo(K))

  # Remove legend
  if (!legend) plt <- plt + ggplot2::theme(legend.position = "none")

  return(plt)
}

#' Helper function for TESS barplots using ggplot
#'
#' @param dat Q matrix
#'
#' @return barplot with Q-values and individuals, colorized by K-value
#'
#' @family TESS functions
#' @export
ggbarplot_helper <- function(gg_df) {
  gg_df %>%
    ggplot2::ggplot(ggplot2::aes(x = order, y = Q_value, fill = cluster)) +
    ggplot2::geom_bar(stat = "identity", col = "darkgray", size = 0.3) +
    ggplot2::scale_y_continuous(expand = c(0,0)) +
    ggplot2::scale_x_discrete(expand = c(0,0)) +
    ggplot2::ylab("Q") +
    ggplot2::labs(fill = "cluster") +
    ggplot2::theme(axis.line = ggplot2::element_line(colour = "black"),
                   axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_blank(),
                   panel.border = ggplot2::element_rect(fill = NA, colour = "black", linetype = "solid", linewidth = 1.5),
                   strip.text.y = ggplot2::element_text(size = 30, face = "bold"),
                   strip.background = ggplot2::element_rect(colour = "white", fill = "white"),
                   panel.spacing = ggplot2::unit(-0.1, "lines"))
}

# Kriging functions -------------------------------------------------------
# -------------------------------------------------------------------------

#' Implements `Krig()` function in the fields package
#'
#' @param qmat matrix of Q-values for given K
#' @param coords matrix of x y coordinates
#' @param grid raster upon which to krige
#'
#' @return
#' @export
new_tess_krig <- function(qmat, coords, grid) {
  # Check CRS
  crs_check(coords, grid)

  # Get K
  K <- ncol(qmat)

  krig_df <- coords_to_sp(coords)

  # Make grid for kriging
  if (!inherits(grid, "SpatRaster")) grid <- terra::rast(grid)
  # krig_grid <- raster_to_grid(grid) # not necessary for type of kriging

  # Krige each K value
  krig_admix <-
    purrr::map(1:K, new_tess_krig_helper, krig_df, qmat, coords, grid) %>%
    terra::rast()

  # mask with original raster layer because the grid fills in all NAs
  # (note: we don't remove NAs because it can change the extent)
  grid <- terra::resample(grid, krig_admix[[1]])
  krig_admix <- terra::mask(krig_admix, grid)

  # Rename layers
  names(krig_admix) <- paste0("K", 1:K)

  # Return stack
  return(krig_admix)
}

#' Krige one K value
#'
#' @param K single K value to krige
#' @param qmat matrix of Q-values
#' @param krig_df df with coords
#' @param coords coordinates in x and y
#' @param grid grid upon which to krige
#'
#' @return
#' @export
new_tess_krig_helper <- function(K, krig_df, qmat, coords, grid) {
  # Add Q values to spatial dataframe
  krig_df$Q <- qmat[, K]

  # Skip if all of the Q values are identical (kriging not possible)
  # if (length(unique(krig_df$Q)) == 1) {
  #   warning(paste0("Only one unique Q value for K = ", K, ", returning NULL (note: may want to consider different K value)"))
  #   return(NULL)
  # }

  co <- capture.output(model <- fields::Krig(coords, qmat[, K]))
  interpolated <- terra::interpolate(raster(grid), model)
  interpolated <- terra::rast(interpolated)
  return(interpolated)
}

# Reproject the coordinates
st_as_sf(coords, coords = c("x", "y"), crs = 4326)
st_transform(coords, 3310)
