format_dist_helper <- function(file_name, output_name){
  path <- here("analysis", "ibdibe", "outputs")
  d <- read.table(here(path, paste0(file_name, ".mdist")))
  id <- read.table(here(path, paste0(file_name, ".mdist.id")))[,2]
  rownames(d) <- colnames(d) <- id
  output_path <- here(path, paste0(output_name, ".csv"))
  write.csv(d, output_path)
  message("wrote dist file to: ", output_path)
}

format_dist <- function(file_name = "58-Sceloporus_annotated_pruned_0.6_chr"){
  message("Formatting ", file_name)
  format_dist_helper(file_name, "58-Sceloporus_dist")
  message("Formatting non-synonymous SNPs")
  format_dist_helper("nonsyn", "nonsyn_dist")
}

get_gendist <- function(){
  path <- here("analysis", "ibdibe", "outputs")
  gendist <- read.csv(here(path, "58-Sceloporus_dist.csv"), row.names = 1)
  colnames(gendist) <- row.names(gendist)
  return(gendist)
}

get_genedist <- function(){
  path <- here("analysis", "ibdibe", "outputs")
  gendist <- read.csv(here(path, "genes_dist.csv"), row.names = 1)
  colnames(gendist) <- row.names(gendist)
  return(gendist)
}

unfold <- function(X, scale = TRUE) {
  x <- vector()
  for (i in 2:nrow(X)) x <- c(x, X[i, 1:i - 1])
  if (scale == TRUE) x <- scale(x, center = TRUE, scale = TRUE)
  return(x)
}

dist_to_df <- function(Y, X, stdz = TRUE){
  #Unfold X and Y
  y <- unfold(Y, scale = stdz)
  dfX <- purrr::map_dfc(X, unfold, scale = stdz) %>% purrr::map_dfc(as.numeric)

  # Make single variable dataframe
  df <- dfX %>%
    dplyr::mutate(Y = y) %>%
    tidyr::gather("var", "X", -Y)

  return(df)
}

gdm_plot_isplines <- function(gdm_model, env, coords, scales = "free", nrow = NULL, ncol = NULL) {
  if (!inherits(coords, "sf")) coords <- sf::st_as_sf(coords, coords = c("x", "y"), crs = sf::st_crs(env))
  splineDat <- gdm::isplineExtract(gdm_model)

  gdm_spline_df <- dplyr::bind_rows(
    data.frame(splineDat$x, var = "x", ID = seq_len(nrow(splineDat$x))),
    data.frame(splineDat$y, var = "y", ID = seq_len(nrow(splineDat$y)))
  ) %>%
    tidyr::pivot_longer(-c("var", "ID"), names_to = "name", values_to = "value") %>%
    tidyr::pivot_wider(names_from = "var", values_from = "value") %>%
    dplyr::mutate(name = factor(name, levels = unique(name))) %>%
    dplyr::mutate(name = gsub("_", " ", name)) %>%
    dplyr::mutate(name = gsub("Geographic", "Geographic distance", name)) %>%
    dplyr::mutate(x = ifelse(name == "Geographic distance", x/1000, x))

  rug_df <- env %>%
    tidyr::pivot_longer(dplyr::everything(), names_to = "name", values_to = "x") %>%
    dplyr::mutate(name = gsub("_", " ", name)) %>%
    dplyr::bind_rows(data.frame(name = "Geographic distance", x = as.vector(sf::st_distance(coords))/1000))

  yend <- min(gdm_spline_df$y, na.rm = TRUE)
  ystart <- yend - (max(gdm_spline_df$y, na.rm = TRUE) * 0.05)

  ggplot2::ggplot(gdm_spline_df) +
    ggplot2::geom_segment(data = rug_df, ggplot2::aes(x = x, xend = x, y = ystart, yend = yend, col = name), alpha = 0.1) +
    ggplot2::geom_line(ggplot2::aes(x = x, y = y, col = name), linewidth = 1) +
    ggplot2::facet_wrap(~name, scales = scales, nrow = nrow, ncol = ncol, strip.position = "bottom") +
    ggplot2::theme_bw() +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::ylab("Partial Regression Distance") +
    ggplot2::xlab("") +
    ggplot2::theme(strip.placement = "outside", strip.background = ggplot2::element_blank(), panel.grid = ggplot2::element_blank(), strip.text = ggplot2::element_text(size = 15), axis.title.y = ggplot2::element_text(size = 14), axis.text = ggplot2::element_text(size = 12), legend.position = "none") +
    ggplot2::scale_color_manual(values = c("Geographic distance" = "black", setNames(scales::hue_pal()(length(setdiff(unique(gdm_spline_df$name), "Geographic distance")) + 1)[-1], setdiff(unique(gdm_spline_df$name), "Geographic distance"))))
}


# convert from matrix/data.frame/sf to formatted df
coords_to_df <- function(coords) {
  if (inherits(coords, "sf")) coords <- sf::st_coordinates(coords)
  if (is.matrix(coords)) coords <- data.frame(coords)
  colnames(coords) <- c("x", "y")
  return(coords)
}

gdm_pc <- function(gdm_model, envlayers, coords, scl = 1, display_axes = FALSE) {
  # convert envlayers to SpatRaster
  if (!inherits(envlayers, "SpatRaster")) envlayers <- terra::rast(envlayers)

  # convert coords to df
  coords <- coords_to_df(coords)

  # CHECK that all of the model variables are included in the stack of environmental layers
  # Create list of environmental predictors (everything but Geographic)
  check_geo <- gdm_model$predictors == "Geographic"
  if (any(check_geo)) {
    model_vars <- gdm_model$predictors[-which(check_geo)]
  } else {
    model_vars <- gdm_model$predictors
  }

  # Check that model variables are included in names of envlayers
  var_check <- model_vars %in% names(envlayers)

  # Print error with missing layers
  if (!all(var_check)) {
    stop(paste("missing model variable(s) from raster stack:", model_vars[!var_check]))
  }

  # Subset envlayers to only include variables in final model
  envlayers_sub <- terra::subset(envlayers, model_vars)

  # CREATE MAP ----------------------------------------------------------------------------------------------------

  # Transform environmental layers
  # TEMPORARY: In new versions of GDM, the input/output rasters are SpatRasters, but for the old version they are rasters
  if (packageVersion("gdm") >= "1.6.0-4") {
    rastTrans <- gdm::gdm.transform(gdm_model, envlayers_sub)
  } else {
    envlayers_sub_raster <- raster::stack(envlayers_sub)
    rastTrans <- gdm::gdm.transform(gdm_model, envlayers_sub_raster)
    rastTrans <- terra::rast(rastTrans)
  }

  # Remove NA values
  rastDat <- na.omit(terra::values(rastTrans))

  # Run PCA
  pcaSamp <- stats::prcomp(rastDat)

  # Count number of layers
  n_layers <- terra::nlyr(rastTrans)
  # Max number of layers to plot is 3, so adjust n_layers accordingly
  if (n_layers > 3) {
    n_layers <- 3
  }

  # Check if there are only coordinate layers
  # If there are only coordinate layers (i.e., no env layers) than you cannot create the map
  if (all(names(rastTrans) %in% c("xCoord", "yCoord"))) stop("All model splines for environmental variables are zero")

  # Make PCA raster
  pcaRast <- terra::predict(rastTrans, pcaSamp, index = 1:n_layers)

  # Scale rasters to get colors (each layer will correspond with R, G, or B in the final plot)
  pcaRastRGB <- stack_to_rgb(pcaRast)

  # If there are fewer than 3 n_layers (e.g., <3 variables), the RGB plot won't work (because there isn't an R, G, and B)
  # To get around this, create a blank raster (i.e., a white raster), and add it to the stack
  if (n_layers < 3) {
    warning("Fewer than three non-zero coefficients provided, adding white substitute layers to RGB plot")
    # Create white raster by multiplying a layer of pcaRast by 0 and adding 255
    white_raster <- pcaRastRGB[[1]] * 0 + 255
  }

  # If n_layers = 2, you end up making a bivariate map
  if (n_layers == 2) {
    pcaRastRGB <- c(pcaRastRGB, white_raster)
  }

  # If n_layers = 1, you end up making a univariate map
  if (n_layers == 1) {
    pcaRastRGB <- c(pcaRastRGB, white_raster, white_raster)
  }

  # Plot variable vectors
  plt <- gdm_pc_helper(pcaSamp, pcaRast, pcaRastRGB, coords, x = "PC1", y = "PC2", scl = scl, display_axes = display_axes)

  if ((n_layers != 3)) {
    stop("variable vector plot is not available for model with fewer than 3 final variables")
  }

  return(plt)
}


gdm_pc_helper <- function(pcaSamp, pcaRast, pcaRastRGB, coords, x = "PC1", y = "PC2", scl = 1, display_axes = FALSE) {
  # Confirm there are exactly 3 axes
  if (terra::nlyr(pcaRastRGB) > 3) {
    stop("Only three PC layers (RGB) can be used for creating the variable plot (too many provided)")
  }
  if (terra::nlyr(pcaRastRGB) < 3) {
    stop("Need exactly three PC layers (RGB) for creating the variable plot (too few provided)")
  }

  # GET PCA DATA ----------------------------------------------------------------------------------------------------

  # Make data frame from PC results
  xpc <- data.frame(pcaSamp$x[, 1:3])

  # Get variable rotations
  varpc <- data.frame(varnames = rownames(pcaSamp$rotation), pcaSamp$rotation)
  varpc <- 
    varpc %>% 
    mutate(varnames = 
      case_when(
        varnames == "xCoord" ~ "X coordinate",
        varnames == "yCoord" ~ "Y coordinate",
        TRUE ~ varnames)
    ) %>%
    # Replace underscores with spaces for plotting
    mutate(varnames = gsub("_", " ", varnames))

  # Get PC values for each coord
  pcavals <- data.frame(terra::extract(pcaRast, coords, ID = FALSE))
  colnames(pcavals) <- colnames(xpc)

  # Rescale var loadings with individual loadings so they fit in the plot nicely
  scldat <- min(
    (max(pcavals[, y], na.rm = TRUE) - min(pcavals[, y], na.rm = TRUE) / (max(varpc[, y], na.rm = TRUE) - min(varpc[, y], na.rm = TRUE))),
    (max(pcavals[, x], na.rm = TRUE) - min(pcavals[, x], na.rm = TRUE) / (max(varpc[, x], na.rm = TRUE) - min(varpc[, x], na.rm = TRUE)))
  )

  # Additionally use a constant scale val (scl) to shrink the final vectors (again for plotting nicely)
  varpc <- data.frame(varpc,
    v1 = scl * scldat * varpc[, x],
    v2 = scl * scldat * varpc[, y]
  )

  # Normalize loadings to [0,1] for RGB
  scale01 <- function(x) {
    rng <- range(x, na.rm = TRUE)
    if (diff(rng) == 0) return(rep(0.5, length(x)))
    (x - rng[1]) / diff(rng)
  }

  varpc <- varpc %>%
    mutate(
      R = scale01(PC1),
      G = scale01(PC2),
      B = scale01(PC3),
      rgb_col = rgb(R, G, B)
    )


  # GET RGB VALS FOR EACH COORD----------------------------------------------------------------------------------------

  pcavalsRGB <- data.frame(terra::extract(pcaRastRGB, coords, ID = FALSE))
  colnames(pcavalsRGB) <- colnames(xpc)

  # Create vector of RGB colors for plotting
  pcacols <- apply(pcavalsRGB, 1, create_rgb_vec)

  # GET RGB VALS FOR ENTIRE RASTER-------------------------------------------------------------------------------------

  # Get sample
  s <- sample(1:terra::ncell(pcaRast), 10000)

  # Get all PC values from raster and remove NAs
  rastvals <- data.frame(terra::values(pcaRast))[s, ]
  colnames(rastvals) <- colnames(xpc)
  rastvals <- rastvals[stats::complete.cases(rastvals), ]

  # Get all RGB values from raster and remove NAs
  rastvalsRGB <- data.frame(terra::values(pcaRastRGB))[s, ]
  colnames(rastvalsRGB) <- colnames(rastvals)
  rastvalsRGB <- rastvalsRGB[stats::complete.cases(rastvalsRGB), ]

  # Create vector of RGB colors for plotting
  rastpcacols <- apply(rastvalsRGB, 1, create_rgb_vec)

  # FINAL PLOT----------------------------------------------------------------------------------------------------------

  # Build base plot
  # Plot points colored by RGB with variable vectors
  plot <- ggplot2::ggplot() +

    # Create axes that cross through origin
    {
      if (display_axes) ggplot2::geom_hline(yintercept = 0, size = 0.2, col = "gray")
    } +
    {
      if (display_axes) ggplot2::geom_vline(xintercept = 0, size = 0.2, col = "gray")
    } +

    # Plot points from entire raster
    #ggplot2::geom_point(data = rastvals, ggplot2::aes_string(x = x, y = y), col = rastpcacols, size = 4, alpha = 0.02) +

    # Plot coord values
    ggplot2::geom_point(data = pcavals, ggplot2::aes_string(x = x, y = y), fill = pcacols, col = "black", pch = 21, size = 2) +

    # Plot variable vectors
    ggplot2::geom_segment(data = varpc, ggplot2::aes(x = 0, y = 0, xend = v1, yend = v2, col = varnames), arrow = ggplot2::arrow(length = ggplot2::unit(0.3, "cm")), linewidth = 0.8, alpha = 1) +
    ggrepel::geom_label_repel(
      data = varpc,
      aes(x = v1, y = v2, label = varnames, col = varnames),
      size = 4,

      # white background behind text
      fill = "white",
      label.size = 0.25,      # border thickness (0 = no border)
      label.r = unit(0.15, "lines"),  # rounded corners

      # breathing room from points
      box.padding = 0.8,
      point.padding = 0.8,

      # Colors
      segment.color = NA
    ) +

    # Plot formatting
    ggplot2::coord_equal() +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.border = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      aspect.ratio = 1, 
      legend.position = "none"
    ) +

  # Add manual color scale so that Latitude and Longitude are always black
  ggplot2::scale_color_manual(
    values = c(
      "X coordinate"  = "black",
      "Y coordinate" = "black",
      setNames(
        scales::hue_pal()(
          length(setdiff(unique(varpc$varnames),
                        c("X coordinate", "Y coordinate"))) + 1
        )[-1],   # ← skip first default ggplot color
        setdiff(unique(varpc$varnames), c("X coordinate", "Y coordinate"))
      )
    )
  )

  # Build plot without PC axes displayed
  if (display_axes == FALSE) {
    plot <- plot +
      # Remove axes
      ggplot2::theme(
        axis.title = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank()
      )
  }

  # Plot
  return(plot)
}

gdm_plot_diss <- function(gdm_model, type = NULL, line = TRUE) {
  obs <- tidyr::as_tibble(gdm_model$observed) %>% dplyr::rename(observed = value)
  pred <- tidyr::as_tibble(gdm_model$predicted) %>% dplyr::rename(predicted = value)
  ecol <- tidyr::as_tibble(gdm_model$ecological) %>% dplyr::rename(ecological = value)
  dat <- dplyr::bind_cols(obs, pred, ecol)
  n <- nrow(dat)

  overlayX_ecol <- seq(min(dat$ecological), max(dat$ecological), length.out = n)
  overlayY_ecol <- 1 - exp(-overlayX_ecol)
  overlayX_pred <- overlayY_pred <- seq(min(dat$predicted), max(dat$predicted), length.out = n)

  plot_ecol <- ggplot2::ggplot(dat) +
    ggplot2::geom_hex(ggplot2::aes(x = ecological, y = observed)) +
    ggplot2::theme_classic() +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::xlab("Predicted ecological distance") +
    ggplot2::ylab("Observed dissimilarity") +
    ggplot2::labs(fill = "Sample\ncount") +
    ggplot2::theme(axis.title = ggplot2::element_text(size = 14), axis.text = ggplot2::element_text(size = 12))

  if (line) plot_ecol <- plot_ecol + ggplot2::geom_line(data = data.frame(x = overlayX_ecol, y = overlayY_ecol), ggplot2::aes(x = x, y = y), color = "tomato2", size = 1)

  plot_pred <- ggplot2::ggplot(dat) +
    ggplot2::geom_hex(ggplot2::aes(x = predicted, y = observed)) +
    ggplot2::theme_classic() +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::xlab("Predicted dissimilarity") +
    ggplot2::ylab("Observed dissimilarity") +
    ggplot2::labs(fill = "Sample\ncount") +
    ggplot2::theme(axis.title = ggplot2::element_text(size = 14), axis.text = ggplot2::element_text(size = 12), legend.title = ggplot2::element_text(size = 14), legend.text = ggplot2::element_text(size = 12))

  if (line) plot_pred <- plot_pred + ggplot2::geom_line(data = data.frame(x = overlayX_pred, y = overlayY_pred), ggplot2::aes(x = x, y = y), color = "tomato2", size = 1)

  legend <- cowplot::get_legend(plot_pred)
  plot_pred <- plot_pred + ggplot2::theme(legend.position = "none")
  plot_ecol <- plot_ecol + ggplot2::theme(legend.position = "none")

  plt <- cowplot::plot_grid(cowplot::plot_grid(plot_ecol, plot_pred, nrow = 1), legend, nrow = 1, rel_widths = c(1, 0.15))

  if (is.null(type)) return(plt)
  if (type == "ecological") return(plot_ecol + ggplot2::theme(legend.position = "right"))
  if (type == "predicted") return(plot_pred + ggplot2::theme(legend.position = "right"))
}
