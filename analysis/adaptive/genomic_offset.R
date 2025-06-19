# Calculating offset functions --------------------------------------------

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
#' @param scale_env attr for scaled env vars
#' @param center_env attr for centered env vars
#'
#' @return list with five elements: projected present, future, offset, global offset, and weights
#' @export
genomic_offset <- function(loadings, biplot, eig, K = 2, env_pres, env_fut) {
  # Extract values from our environmental rasters
  env <- terra::extract(env_pres, coords)
  # Standardize environmental variables and make into dataframe
  env <- scale(env, center = TRUE, scale = TRUE)
  # Recovering scaling coefficients for extracted env values
  scale_env <- attr(env, 'scaled:scale')
  center_env <- attr(env, 'scaled:center')

  # Deal with future layer naming; 1 is RCP2.6 and 2 is RCP8.5
  env_fut_26 <- terra::subset(env_fut, c(1,3)) # BIO1 ssp126 & NDVI
  names(env_fut_26) <- names(env_pres)
  env_fut_85 <- terra::subset(env_fut, 2:3) # BIO ssp585 & NDVI
  names(env_fut_85) <- names(env_pres)

  var_env_proj_pres <- offset_scaling_helper(env_layer = env_pres, center_env, scale_env, biplot)
  var_env_proj_fut_26 <- offset_scaling_helper(env_layer = env_fut_26, center_env, scale_env, biplot)
  var_env_proj_fut_85 <- offset_scaling_helper(env_layer = env_fut_85, center_env, scale_env, biplot)

  Proj_pres <- offset_proj_helper(biplot, var_env_proj = var_env_proj_pres, K, type = "present")
  Proj_fut_26 <- offset_proj_helper(biplot, var_env_proj = var_env_proj_fut_26, K, type = "future")
  Proj_fut_85 <- offset_proj_helper(biplot, var_env_proj = var_env_proj_fut_85, K, type = "future")

pdf(paste0(here("analysis", "adaptive", "plots"), "/test.pdf"), width = 12, height = 8)
plot(Proj_offset_pres[[1]])
dev.off()


  
  # Single axis genetic offset
  Proj_offset_26 <- list()
  for(i in 1:K) {
    Proj_offset_26[[i]] <- abs(Proj_pres[[i]] - Proj_fut_26[[i]])
    names(Proj_offset_26)[i] <- paste0("RDA", as.character(i))
  }
  Proj_offset_85 <- list()
  for(i in 1:K) {
    Proj_offset_85[[i]] <- abs(Proj_pres[[i]] - Proj_fut_85[[i]])
    names(Proj_offset_85)[i] <- paste0("RDA", as.character(i))
  }

  # Weights based on axis eigenvalues
  weights <- eig %>% dplyr::mutate(weights = mod.CCA.eig / sum(eig$mod.CCA.eig)) %>% pull(weights)

  # Weighing the current and future adaptive indices based on the eigenvalues of the associated axes
  # Proj_offset_pres <- do.call(cbind, lapply(1:K, function(x) rasterToPoints(Proj_pres[[x]])[,-c(1,2)]))
  Proj_offset_pres <- do.call(cbind, lapply(1:K, function(x) {terra::values(Proj_pres[[x]], mat = TRUE)}))
  Proj_offset_pres <- as.data.frame(do.call(cbind, lapply(1:K, function(x) Proj_offset_pres[,x] * weights[x])))

  # Proj_offset_fut_26 <- do.call(cbind, lapply(1:K, function(x) rasterToPoints(Proj_fut_26[[x]])[,-c(1,2)]))
  Proj_offset_fut_26 <- do.call(cbind, lapply(1:K, function(x) {terra::values(Proj_fut_26[[x]], mat = TRUE)}))
  Proj_offset_fut_26 <- as.data.frame(do.call(cbind, lapply(1:K, function(x) Proj_offset_fut_26[,x] * weights[x])))

  # Proj_offset_fut_85 <- do.call(cbind, lapply(1:K, function(x) rasterToPoints(Proj_fut_2[[x]])[,-c(1,2)]))
  Proj_offset_fut_85 <- do.call(cbind, lapply(1:K, function(x) {terra::values(Proj_fut_85[[x]], mat = TRUE)}))
  Proj_offset_fut_85 <- as.data.frame(do.call(cbind, lapply(1:K, function(x) Proj_offset_fut_85[,x] * weights[x])))

  # Predict a global genetic offset, incorporating the first K axes weighted by their eigenvalues
  # TODO in development
  # ras_26 <- Proj_offset_26[[1]]
  # ras_26[!is.na(ras_26)] <- unlist(lapply(1:nrow(Proj_offset_pres), function(x) dist(rbind(Proj_offset_pres[x,], Proj_offset_fut_26[x,]), method = "euclidean")))
  # names(ras_26) <- "Global_offset_26"
  # Proj_offset_global_26 <- ras_26

  # ras_85 <- Proj_offset_85[[1]]
  # ras_85[!is.na(ras_85)] <- unlist(lapply(1:nrow(Proj_offset_pres), function(x) dist(rbind(Proj_offset_pres[x,], Proj_offset_fut_85[x,]), method = "euclidean")))
  # names(ras_85) <- "Global_offset_85"
  # Proj_offset_global_85 <- ras_85

  # Return projections for current and future climates for each RDA axis, prediction of genetic offset for each RDA axis and a global genetic offset
  return(list(Proj_pres = Proj_pres, Proj_fut_RCP26 = Proj_fut_26, Proj_fut_RCP85 = Proj_fut_85,
              Proj_offset_RCP26 = Proj_offset_26, Proj_offset_RCP85 = Proj_offset_85,
              Proj_offset_global_RCP26 = NULL, Proj_offset_global_RCP85 = NULL,
              weights = weights[1:K]))
}

offset_scaling_helper <- function(env_layer, center_env, scale_env, biplot) {
  # Matrix with NA values preserved
  env_vals <- terra::values(env_layer[[row.names(biplot)]], mat = TRUE)
  # Manually scale the matrix using precomputed center and scale
  scaled_vals <- sweep(env_vals, 2, center_env[row.names(biplot)], "-")
  scaled_vals <- sweep(scaled_vals, 2, scale_env[row.names(biplot)], "/")
  # Convert to df
  var_env_proj <- as.data.frame(scaled_vals)
  return(var_env_proj)
}

#' Helper function for calculating offset
#'
#' @param biplot RDA biplot results
#' @param var_env_proj projected env
#' @param K number of layers
#' @param type either "present" or "future"; just for naming
#'
#' @return
#' @export
offset_proj_helper <- function(biplot, var_env_proj, K, type) {
  Proj_list <- list()
  if (type == "present") prefix = "RDA_pres_"
  if (type == "future") prefix = "RDA_fut_"
  for(i in 1:K) {
    vars_i <- rownames(biplot)
    loadings_i <- biplot[, i]
  
    projection_vals <- rowSums(var_env_proj[, vars_i] * matrix(loadings_i, nrow = nrow(var_env_proj), ncol = length(loadings_i), byrow = TRUE))
    ras <- env_pres[[1]]
    ras[] <- projection_vals

    names(ras) <- paste0(prefix, as.character(i))
    Proj_list[[i]] <- ras
    names(Proj_list)[i] <- paste0("RDA", as.character(i))
  }

  # Old, uses raster package
  # for(i in 1:K) {
  #   ras <- env[[1]]
  #   ras[!is.na(ras)] <- as.vector(apply(var_env_proj[,rownames(biplot[i])], 1, function(x) sum(x * biplot[,i])))
  #   names(ras) <- paste0(prefix, as.character(i))
  #   Proj_list[[i]] <- ras
  #   names(Proj_list)[i] <- paste0("RDA", as.character(i))
  # }
  return(Proj_list)
}

#' @title Predict Biological Dissimilarities Between Sites or Times Using a
#' Fitted Generalized Dissimilarity Model
#'
#' @description This function predicts biological distances between sites or times using a
#'  model object returned from \code{\link[gdm]{gdm}}. Predictions between site
#'  pairs require a data frame containing the values of predictors for pairs
#'  of locations, formatted as follows: distance, weights, s1.X, s1.Y, s2.X,
#'  s2.Y, s1.Pred1, s1.Pred2, ..., s1.PredN, s2.Pred1, s2.Pred2, ..., s2.PredN, ...,
#'  Predictions of biological change through time require two raster stacks or
#'  bricks for environmental conditions at two time periods, each with a
#'  layer for each environmental predictor in the fitted model.
#'  
#'  https://github.com/fitzLab-AL/gdm/blob/master/R/gdm.predict.R
#'
#' @usage \method{predict}{gdm}(object, data, time=FALSE, predRasts=NULL, filename="", ...)
#'
#' @param object A gdm model object resulting from a call to \code{\link[gdm]{gdm}}.
#'
#' @param data Either a data frame containing the values of predictors for pairs
#' of sites, in the same format and structure as used to fit the model using
#' \code{\link[gdm]{gdm}}, or a raster stack if a prediction of biological change
#' through time is needed.
#'
#' For a data frame, the first two columns - distance and weights - are required
#' by the function but are not used in the prediction and can therefore be filled
#' with dummy data (e.g. all zeros). If geo is TRUE, then the s1.X, s1.Y and s2.X,
#' s2.Y columns will be used for calculating the geographical distance between
#' each site for inclusion of the geographic predictor term into the GDM model.
#' If geo is FALSE, then the s1.X, s1.Y, s2.X and s2.Y data columns are ignored.
#' However these columns are still REQUIRED and can be filled with dummy data
#' (e.g. all zeroes). The remaining columns are for N predictors for Site 1 and
#' followed by N predictors for Site 2. The order of the columns must match those
#' in the site-pair table used to fit the model.
#'
#' A raster stack should be provided only when time=T and should contain one
#' layer for each environmental predictor in the same order as the columns in
#' the site-pair table used to fit the model.
#'
#' @param time TRUE/FALSE: Is the model prediction for biological change through time?
#' @param predRasts A raster stack characterizing environmental conditions for a
#' different time in the past or future, with the same extent, resolution, and
#' layer order as the data object. Required only if time=T.
#' @param filename character. Output filename for rasters. When provided the raster layers are
#' written to file directly.
#' @param ... additional arguments to pass to terra \code{\link[terra]{predict}} function.
#'
#' @return predict returns either a response vector with the same length as the
#'  number of rows in the input data frame or a raster depicting change through time across the study region.
#'
predict.gdm <- function(object, data, time = FALSE, predRasts = NULL, filename = "", ...) {
  # object <- object$model
  if (time) {
    for(i in 1:terra::nlyr(data)){
      if(names(data)[i] != names(predRasts)[i]){
        stop("Layer names do not match the variables used to fit the model.")
      }
    }
    if(terra::nlyr(data) != length(object$predictors) - 1 | terra::nlyr(predRasts) != length(object$predictors) - 1){
        stop("Number of variables supplied for prediction does not equal the number used to fit the model.")
      }
    # if (object$geo) {
    #   if(terra::nlyr(data) != length(object$predictors) - 1 | terra::nlyr(predRasts) != length(object$predictors) - 1){
    #     stop("Number of variables supplied for prediction does not equal the number used to fit the model.")
    #   }
    # } else {
    #   if(terra::nlyr(data)!=length(object$predictors) | terra::nlyr(predRasts)!=length(object$predictors)){
    #     stop("Number of variables supplied for prediction does not equal the number used to fit the model.")
    #   }
    # }

    # create XY rasters; data and predRasts must have the same XY
    x <- terra::init(data[[1]], fun = "x")
    y <- terra::init(data[[1]], fun = "y")
    
    # sets the correct names to the data
    names(data) <- paste0("s1.", names(data))
    names(predRasts) <- paste0("s2.", names(predRasts))
    
    # stack all the raster layers to for prediction
    data <- c(
      # stats::setNames(dummData, "distance"),
      # stats::setNames(dummData, "weights"),
      stats::setNames(x, "s1.xCoord"),
      stats::setNames(y, "s1.yCoord"),
      stats::setNames(x, "s2.xCoord"),
      stats::setNames(y, "s2.yCoord"),
      data,
      predRasts
    )
    
  }
  
  # makes the prediction based on the data object
  gdm_predict <- function(mod, dat, ...) {
    nr <- nrow(dat)
    predicted <- rep(0, times = nr)
    
    # convert to matrix once
    dat <- as.matrix(dat)
    # add the constants
    const <- matrix(0L, nrow = nr, ncol = 2)
    colnames(const) <- c("distance", "weights")
    dat <- cbind(const, dat)
    
    z <- .C( "GDM_PredictFromTable",
             dat,
             as.integer(mod$geo),
             as.integer(length(mod$predictors)),
             as.integer(nr),
             as.double(mod$knots),
             as.integer(mod$splines),
             as.double(c(mod$intercept, mod$coefficients)),
             preddata = as.double(predicted),
             PACKAGE = "gdm")
    
    return(z$preddata)
  }
  
  # if a time prediction, maps the predicted values to a raster and returns
  # the layer, otherwise returns a dataframe of the predicted values
  if (time) {
    # predict using gdm model and terra package
    output <- terra::predict(
      object = data,
      model = object,
      fun = gdm_predict,
      na.rm = TRUE,
      filename = here("analysis", "adaptive", "outputs", "58-Sceloporus_GDMoffset.tif"), overwrite = TRUE
    )
    
    return(output)
    
  } else {
    # predict using a data.frame
    output <- gdm_predict(
      mod = object,
      dat = data
    )
    
    return(output)
  }
}

# Plotting offset functions -----------------------------------------------

#' Plot genomic offset maps
#'
#' @param env offset rasters
#' @param plot_type options are "basic", "rainbow", "extracted_vals", or "extracted_rainbow"
#' @param bkg shape file for background plotting
#' @param free_scales whether to have separate scales or not for RDA axes
#' @param index_name name for legend
#' @param viridis_option option for viridis coloring (defaults to "B")
#' @param coords if `plot_type = "extracted_vals" or "extracted_rainbow"`, sampling coordinates
#' @param bkg_col if `plot_type = "rainbow"`, background color for blank raster if fewer than 3 layers
#'
#' @returns
#' @export
plot_offset <- function(env, bkg, plot_type = "basic", free_scales = FALSE,
                        index_name = "Genomic offset", viridis_option = "B", coords = NULL,
                        biplot_axes = c(1, 2), bkg_col = "white") {
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
                                  loadings = loadings, biplot_axes = biplot_axes, coords = coords, bkg_col = bkg_col)
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
rainbow_map_offset <- function(Proj_data, bkg, n_layers, loadings, biplot_axes, coords, bkg_col) {
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
    if (bkg_col == "white") bkg_raster <- aiRGB[[1]] * 0 + 255
    if (bkg_col == "black") bkg_raster <- aiRGB[[1]] * 0
  }

  # If n_layers = 2, you end up making a bivariate map
  if (n_layers == 2) {
    aiRGB <- c(aiRGB[[1]], aiRGB[[2]], bkg_raster)
  }

  # If n_layers = 1, you end up making a univariate map
  if (n_layers == 1) {
    aiRGB <- c(aiRGB, bkg_raster, bkg_raster)
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
  plot <-
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
