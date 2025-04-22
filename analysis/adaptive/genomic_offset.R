#' Predict genomic offset from an RDA model
#' Code adapted from Capblancq & Forester (2021) https://doi.org/10.1111/2041-210X.13722
#' GitHub repo available here: https://github.com/Capblancq/RDA-landscape-genomics/blob/main/src/genomic_offset.R
#'
#' @param loadings loadings from RDA model
#' @param biplot biplot values from RDA model
#' @param eig eigenvalues from RDA model
#' @param K number of RDA axes to retain (defaults to 2)
#' @param env_pres present env layers
#' @param env_fut future env layers
#' @param range species range; if provided offset predictions will be masked to range
#' @param method
#' @param scale_env whether to scale env vars
#' @param center_env whether to center env vars
#' @param mod RDA model; only if `method = "predict"`
#'
#' @return list with five elements: projected present, future, offset, global offset, and weights
#' @export
genomic_offset <- function(loadings, biplot, eig, K = 2, env_pres, env_fut, range = NULL, method = "loadings", scale_env, center_env, mod = NULL) {
  # Mask with the range if supplied
  if(!is.null(range)){
    env_pres <- raster::mask(env_pres, range)
    env_fut <- raster::mask(env_fut, range)
  }

  # Deal with future layer naming; 1 is RCP2.6 and 2 is RCP8.5
  env_fut_1 <- raster::subset(env_fut, c(1,3)) # BIO1 ssp125 & NDVI
  names(env_fut_1) <- names(env_pres)
  env_fut_2 <- raster::subset(env_fut, 2:3) # BIO ssp585 & NDVI
  names(env_fut_2) <- names(env_pres)

  # Formatting and scaling environmental rasters for projection
  var_env_proj_pres <- as.data.frame(scale(raster::rasterToPoints(env_pres[[row.names(biplot)]])[,-c(1,2)], center_env[row.names(biplot)], scale_env[row.names(biplot)]))
  var_env_proj_fut_1 <- as.data.frame(scale(raster::rasterToPoints(env_fut_1[[row.names(biplot)]])[,-c(1,2)], center_env[row.names(biplot)], scale_env[row.names(biplot)]))
  var_env_proj_fut_2 <- as.data.frame(scale(raster::rasterToPoints(env_fut_2[[row.names(biplot)]])[,-c(1,2)], center_env[row.names(biplot)], scale_env[row.names(biplot)]))

  # Predicting pixels genetic component based on the loadings of the variables
  if(method == "loadings"){
    # Projection for each RDA axis; TODO warnings generated?
    Proj_pres <- offset_proj_helper(env = env_pres, biplot = biplot, var_env_proj = var_env_proj_pres, K = K, type = "present")
    Proj_fut_1 <- offset_proj_helper(env = env_fut_1, biplot = biplot, var_env_proj = var_env_proj_fut_1, K = K, type = "future")
    Proj_fut_2 <- offset_proj_helper(env = env_fut_2, biplot = biplot, var_env_proj = var_env_proj_fut_2, K = K, type = "future")

    # Single axis genetic offset
    Proj_offset_1 <- list()
    for(i in 1:K) {
      Proj_offset_1[[i]] <- abs(Proj_pres[[i]] - Proj_fut_1[[i]])
      names(Proj_offset_1)[i] <- paste0("RDA", as.character(i))
    }
    Proj_offset_2 <- list()
    for(i in 1:K) {
      Proj_offset_2[[i]] <- abs(Proj_pres[[i]] - Proj_fut_2[[i]])
      names(Proj_offset_2)[i] <- paste0("RDA", as.character(i))
    }
  }

  # Weights based on axis eigenvalues
  weights <- eig %>% dplyr::mutate(weights = mod.CCA.eig / sum(eig$mod.CCA.eig)) %>% pull(weights)

  # Weighing the current and future adaptive indices based on the eigenvalues of the associated axes
  Proj_offset_pres <- do.call(cbind, lapply(1:K, function(x) rasterToPoints(Proj_pres[[x]])[,-c(1,2)]))
  Proj_offset_pres <- as.data.frame(do.call(cbind, lapply(1:K, function(x) Proj_offset_pres[,x] * weights[x])))

  Proj_offset_fut_1 <- do.call(cbind, lapply(1:K, function(x) rasterToPoints(Proj_fut_1[[x]])[,-c(1,2)]))
  Proj_offset_fut_1 <- as.data.frame(do.call(cbind, lapply(1:K, function(x) Proj_offset_fut_1[,x] * weights[x])))

  Proj_offset_fut_2 <- do.call(cbind, lapply(1:K, function(x) rasterToPoints(Proj_fut_2[[x]])[,-c(1,2)]))
  Proj_offset_fut_2 <- as.data.frame(do.call(cbind, lapply(1:K, function(x) Proj_offset_fut_2[,x] * weights[x])))

  # Predict a global genetic offset, incorporating the first K axes weighted by their eigenvalues
  ras_1 <- Proj_offset_1[[1]]
  ras_1[!is.na(ras_1)] <- unlist(lapply(1:nrow(Proj_offset_pres), function(x) dist(rbind(Proj_offset_pres[x,], Proj_offset_fut_1[x,]), method = "euclidean")))
  names(ras_1) <- "Global_offset_1"
  Proj_offset_global_1 <- ras_1

  ras_2 <- Proj_offset_2[[1]]
  ras_2[!is.na(ras_2)] <- unlist(lapply(1:nrow(Proj_offset_pres), function(x) dist(rbind(Proj_offset_pres[x,], Proj_offset_fut_2[x,]), method = "euclidean")))
  names(ras_2) <- "Global_offset_2"
  Proj_offset_global_2 <- ras_2

  # Return projections for current and future climates for each RDA axis, prediction of genetic offset for each RDA axis and a global genetic offset
  return(list(Proj_pres = Proj_pres, Proj_fut_RCP26 = Proj_fut_1, Proj_fut_RCP85 = Proj_fut_2,
              Proj_offset_RCP26 = Proj_offset_1, Proj_offset_RCP85 = Proj_offset_2,
              Proj_offset_global_RCP26 = Proj_offset_global_1, Proj_offset_global_RCP85 = Proj_offset_global_2,
              weights = weights[1:K]))
}

#' Helper function for calculating offset
#'
#' @param env env layers to project
#' @param biplot RDA biplot results
#' @param var_env_proj projected env
#' @param K number of layers
#' @param type either "present" or "future"; just for naming
#'
#' @return
#' @export
offset_proj_helper <- function(env, biplot, var_env_proj, K, type) {
  Proj_list <- list()
  if (type == "present") prefix = "RDA_pres_"
  if (type == "future") prefix = "RDA_fut_"
  for(i in 1:K) {
    ras <- env[[1]]
    ras[!is.na(ras)] <- as.vector(apply(var_env_proj[,rownames(biplot[i])], 1, function(x) sum(x * biplot[,i])))
    names(ras) <- paste0(prefix, as.character(i))
    Proj_list[[i]] <- ras
    names(Proj_list)[i] <- paste0("RDA", as.character(i))
  }
  return(Proj_list)
}

#' Plot genomic offset maps
#'
#' @param env offset rasters
#' @param plot_type options are "basic", "rainbow", "extracted_vals", or "extracted_rainbow"
#' @param bkg shape file for background plotting
#' @param free_scales whether to have separate scales or not for RDA axes
#' @param index_name name for legend
#' @param viridis_option option for viridis coloring (defaults to "B")
#' @param coords if `plot_type = "extractd_vals" or "extracted_rainbow"`, sampling coordinates
#'
#' @returns
#' @export
plot_offset <- function(env, bkg, plot_type = "basic", free_scales = FALSE,
                        index_name = "Genomic offset", viridis_option = "B", coords = NULL,
                        biplot_axes = c(1, 2)) {
  if (inherits(env, "RasterStack")) env <- terra::rast(env)
  n_layers = nlyr(env)

  if (plot_type == "basic") {
    if (!free_scales) {
      p <- ggplot() +
        geom_sf(data = bkg, fill = "lightgrey") +
        geom_spatraster(data = env) +
        scale_fill_viridis_c(option = viridis_option, na.value = "transparent") +
        geom_sf(data = bkg, fill = NA, size = 0.1) +
        xlab("Longitude") +
        ylab("Latitude") +
        ggplot2::guides(fill = guide_legend(title = paste0(index_name))) +
        facet_grid(~lyr) +
        theme_map() +
        theme(panel.grid = element_blank(), plot.background = element_blank(),
              panel.background = element_blank(), strip.text = element_text(size = 11))
    }
    if (free_scales) {
      p1 <- ggplot() +
        geom_sf(data = bkg, fill = "lightgrey") +
        geom_spatraster(data = env[[1]]) +
        scale_fill_viridis_c(option = viridis_option, na.value = "transparent") +
        geom_sf(data = bkg, fill = NA, size = 0.1) +
        xlab("Longitude") +
        ylab("Latitude") +
        ggplot2::guides(fill = guide_legend(title = paste0(index_name))) +
        theme_map() +
        theme(panel.grid = element_blank(), plot.background = element_blank(),
              panel.background = element_blank(), strip.text = element_text(size = 11))
      p2 <- ggplot() +
        geom_sf(data = bkg, fill = "lightgrey") +
        geom_spatraster(data = env[[2]]) +
        scale_fill_viridis_c(option = viridis_option, na.value = "transparent") +
        geom_sf(data = bkg, fill = NA, size = 0.1) +
        xlab("Longitude") +
        ylab("Latitude") +
        ggplot2::guides(fill = guide_legend(title = paste0(index_name))) +
        theme_map() +
        theme(panel.grid = element_blank(), plot.background = element_blank(),
              panel.background = element_blank(), strip.text = element_text(size = 11))
      if (n_layers == 2) p <- plot_grid(p1, p2, nrow = 1)
      if (n_layers == 3) {
        p3 <- ggplot() +
          geom_sf(data = bkg, fill = "lightgrey") +
          geom_spatraster(data = env[[3]]) +
          scale_fill_viridis_c(option = viridis_option, na.value = "transparent") +
          geom_sf(data = bkg, fill = NA, size = 0.1) +
          xlab("Longitude") +
          ylab("Latitude") +
          ggplot2::guides(fill = guide_legend(title = paste0(index_name))) +
          theme_map() +
          theme(panel.grid = element_blank(), plot.background = element_blank(),
                panel.background = element_blank(), strip.text = element_text(size = 11))
        p <- plot_grid(p1, p2, p3, nrow = 1)
      }
    }
  }

  if (plot_type == "extracted_vals") {
    ext <- terra::extract(env, coords, ID = FALSE, xy = TRUE)
    tidy_ext <- ext %>% pivot_longer(cols = 1:n_layers, names_to = "axis", values_to = "value")
    axis_names <- unique(tidy_ext$axis)
    p1 <- ggplot() +
      geom_sf(data = bkg, color = "lightgrey") +
      geom_point(data = tidy_ext %>% filter(axis == axis_names[[1]]),
                 aes(x = x, y = y, color = value), size = 3) +
      scale_color_viridis_c(option = viridis_option, name = index_name) +
      theme_map()
    p2 <- ggplot() +
      geom_sf(data = bkg, color = "lightgrey") +
      geom_point(data = tidy_ext %>% filter(axis == axis_names[[2]]),
                 aes(x = x, y = y, color = value), size = 3) +
      scale_color_viridis_c(option = viridis_option, name = index_name) +
      theme_map()

    if (n_layers == 2) p <- plot_grid(p1, p2, nrow = 1)
    if (n_layers == 3) {
      p3 <- ggplot() +
        geom_sf(data = bkg, color = "lightgrey") +
        geom_point(data = tidy_ext %>% filter(axis == axis_names[[3]]),
                   aes(x = x, y = y, color = value), size = 3) +
        scale_color_viridis_c(option = viridis_option, name = index_name) +
        theme_map()
      p <- plot_grid(p1, p2, p3, nrow = 1)
    }
  }

  if (plot_type == "rainbow" | plot_type == "extracted_rainbow") {
    rainbow <- rainbow_map_offset(Proj_data = env, bkg = bkg, n_layers = n_layers,
                                  loadings = loadings, biplot_axes = biplot_axes, coords = coords)
    if (plot_type == "rainbow") p <- plot_grid(rainbow$map, rainbow$vector_load, rel_widths = c(2, 1))
    if (plot_type == "extracted_rainbow") {
      p1 <- ggplot() +
        geom_sf(data = bkg, color = "lightgrey") +
        geom_point(data = coords, aes(x = x, y = y), fill = rainbow$pcacols, color = "black", pch = 21, size = 3) +
        theme_map()
      p <- plot_grid(p1, rainbow$vector_load, nrow = 1, rel_widths = c(2, 1))
    }
  }
  return(p)
}

#' Build rainbow map with vector loadings for genomic offset
#'
#' @param Proj_data
#' @param bkg
#' @param n_layers
#' @param loadings
#' @param biplot_axes
#' @param coords
#'
#' @returns
#' @export
rainbow_map_offset <- function(Proj_data, bkg, n_layers, loadings, biplot_axes, coords) {
  if (inherits(Proj_data, "RasterStack")) Proj_data <- terra::rast(Proj_data)
  # Max number of layers to plot is 3, so adjust n_layers accordingly
  if (n_layers > 3) {
    n_layers <- 3
  }
  # Scale rasters to get colors (each layer will correspond with R, G, or B in the final plot)
  aiRGB <- stack_to_rgb(Proj_data)

  # If there are fewer than 3 n_layers (e.g., <3 variables), the RGB plot won't work (because there isn't an R, G, and B)
  # To get around this, create a blank raster (i.e., a white raster), and add it to the stack
  if (n_layers < 3) {
    warning("Fewer than three non-zero coefficients provided, adding white substitute layers to RGB plot")
    # Create white raster by multiplying a layer of pcaRast by 0 and adding 255
    white_raster <- aiRGB[[1]] * 0 + 255
  }

  # If n_layers = 2, you end up making a bivariate map
  if (n_layers == 2) {
    aiRGB <- c(aiRGB[[1]], aiRGB[[2]], white_raster)
  }

  # If n_layers = 1, you end up making a univariate map
  if (n_layers == 1) {
    aiRGB <- c(aiRGB, white_raster, white_raster)
  }
  p_rainbow <- ggplot() +
    geom_sf(data = bkg, fill = "lightgrey") +
    geom_spatraster_rgb(data = aiRGB, r = 3, g = 1, b = 2) +
    theme_map()
  p_var <- plot_var_loadings_offset(Proj_data = Proj_data, loadings = loadings,
                             biplot_axes = biplot_axes, aiRGB = aiRGB, coords = coords)
  return(list(map = p_rainbow, vector_load = p_var$var_load_plot, pcacols = p_var$pcacols, rastRGB = aiRGB))
}

#' Plot variable vector loadings as legend for genomic offset rainbow map
#'
#' @param Proj_data
#' @param loadings
#' @param biplot_axes
#' @param aiRGB
#' @param coords
#'
#' @returns
#' @export
plot_var_loadings_offset <- function(Proj_data, loadings, biplot_axes, aiRGB, coords) {
  # TAB_inds <- data.frame(names = rownames(loadings %>% filter(score == "sites")), loadings %>% filter(score == "sites"))
  TAB_var <- data.frame(names = rownames(loadings %>% filter(score == "biplot")), loadings %>% filter(score == "biplot"))

  # Select axes for plotting
  xax <- paste0("RDA", biplot_axes[1], "_offset")
  yax <- paste0("RDA", biplot_axes[2], "_offset")
  # Extract offset values for points
  ext <- terra::extract(Proj_data, coords, ID = FALSE, xy = TRUE)
  if (n_layers == 2) ext <- ext %>% rename("RDA1_offset" = 1, "RDA2_offset" = 2)
  if (n_layers == 3) ext <- ext %>% rename("RDA1_offset" = 1, "RDA2_offset" = 2, "RDA3_offset" = 3)
  ext_sub <- ext %>% dplyr::select(any_of(c(xax, yax)))
  colnames(ext_sub) <- c("x", "y")
  TAB_var_sub <- TAB_var %>% dplyr::select(c("RDA1", "RDA2"))
  colnames(TAB_var_sub) <- c("x", "y")

  # Scale the variable loadings for the arrows
  TAB_var_sub$x <- TAB_var_sub$x * max(ext_sub$x) / stats::quantile(TAB_var_sub$x)[4]
  TAB_var_sub$y <- TAB_var_sub$y * max(ext_sub$y) / stats::quantile(TAB_var_sub$y)[4]

  # GET RGB VALS FOR SAMPLES-------------------------------------------------------------------------------------

  # Get the colors from the rainbow plot for each ind
  pts <- data.frame(terra::extract(aiRGB, coords, ID = FALSE))
  # colnames(pts) <- colnames(xpc)
  # Create vector of RGB colors for plotting
  pcacols <- apply(pts, 1, create_rgb_vec)

  # GET RGB VALS FOR ENTIRE RASTER-------------------------------------------------------------------------------------

  # Get sample
  s <- sample(1:terra::ncell(Proj_data), 10000)

  # Get all PC values from raster and remove NAs
  rastvals <- data.frame(terra::values(Proj_data))[s, ]
  colnames(rastvals) <- c("x", "y")
  rastvals <- rastvals[stats::complete.cases(rastvals), ]

  # Get all RGB values from raster and remove NAs
  rastvalsRGB <- data.frame(terra::values(aiRGB))[s, ]
  colnames(rastvalsRGB) <- colnames(rastvals)
  rastvalsRGB <- rastvalsRGB[stats::complete.cases(rastvalsRGB), ]

  # Create vector of RGB colors for plotting
  rastpcacols <- apply(rastvalsRGB, 1, create_rgb_vec)

  # FINAL PLOT ----------------------------------------
  # Plot points colored by RGB with variable vectors
  ggplot() +
    # geom_hline(yintercept = 0, linewidth = 0.5, col = "gray") +
    # geom_vline(xintercept = 0, linewidth = 0.5, col = "gray") +
    geom_point(data = rastvals, aes(x = x, y = y), col = rastpcacols, size = 4, alpha = 0.02) +
    geom_point(data = ext_sub, aes(x = x, y = y), fill = pcacols, col = "black", pch = 21, size = 3) +
    # geom_segment(data = TAB_var_sub, aes(xend = x, yend = y, x = 0, y = 0), color = "black", linewidth = 0.15, linetype = 1, arrow = ggplot2::arrow(length = ggplot2::unit(0.02, "npc"))) +
    # ggrepel::geom_text_repel(data = TAB_var_sub, aes(x = x, y = y, label = rownames(TAB_var_sub)), size = 4) +
    xlab(xax) +
    ylab(yax) +
    # Plot formatting
    ggplot2::coord_equal() +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.border = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      aspect.ratio = 1
    )
  return(list(var_load_plot = plot, pcacols = pcacols))
}

#' Helper function to create rgb vector
#'
#' @export
#' @noRd
create_rgb_vec <- function(vec) {
  if (any(is.na(vec))) x <- NA else x <- rgb(vec[3], vec[1], vec[2], maxColorValue = 255)
  return(x)
}
