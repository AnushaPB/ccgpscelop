# Import/export functions -------------------------------------------------

#' Import environmental layers of choice for RDA, adaptive index, or genomic offset calculations
#'
#' @param type options are "rasterPCs" or "ind_layers" for BIO1 + NDVI
#' @param future whether to also import future env layers (defaults to FALSE)
#'
#' @return
#' @export
import_env_files <- function(type = "rasterPCs", future = FALSE) {
  if (type == "rasterPCs") {
    env_pres <- terra::rast(here("data", "env", "california_chelsa_bioclim_1981-2010_V.2.1_pca.tif"))
    names(env_pres) <- paste("env_", names(env_pres), sep = "")
  }
  if (type == "ind_layers") {
    all_pres <- terra::rast(here("data", "env", "california_chelsa_bioclim_1981-2010_V.2.1.tif"))
    bio1 <- terra::subset(all_pres, "CHELSA_bio1_1981-2010_V.2.1")
    ndvi <- terra::rast(here("data", "env", "california_ndvi_mean_2000_2020.tif"))
    # Stack layers together and rename
    resamp_ndvi <- terra::resample(ndvi, bio1)
    env_pres <- c(bio1, resamp_ndvi)
    names(env_pres) <- c("BIO1", "NDVI")
  }

  if (future) {
    cap_model <- "GFDL-ESM4" # "IPSL-CM6A-LR"
    lowerc_model <- "gfdl-esm4"
    RCP = c(2.6, 8.5)
    ssp = c("ssp126", "ssp585")

    if (type == "rasterPCs") {
      env_fut_1 <- terra::rast(paste0(here("data", "env", "future"), "/CHELSA_2071-2100_", cap_model, "_", ssp[1], "_V.2.1_pca.tif"))
      env_fut_2 <- terra::rast(paste0(here("data", "env", "future"), "/CHELSA_2071-2100_", cap_model, "_", ssp[2], "_V.2.1_pca.tif"))
      env_fut <- c(env_fut_1, env_fut_2)
    }

    if (type == "ind_layers") {
      bio1_fut_1 <- terra::rast(paste0(here("data", "env", "future", "envicloud/chelsa/chelsa_V2/GLOBAL/climatologies"), "/2071-2100/", cap_model, "/", ssp[1], "/bio/", "CHELSA_bio1_2071-2100_", lowerc_model, "_", ssp[1], "_V.2.1.tif"))
      cropped_1 <- terra::crop(bio1_fut_1, env_pres[[1]], mask = TRUE) # is this necessary?
      resamp_1 <- terra::resample(cropped_1, env_pres[[1]])

      bio1_fut_2 <- terra::rast(paste0(here("data", "env", "future", "envicloud/chelsa/chelsa_V2/GLOBAL/climatologies"), "/2071-2100/", cap_model, "/", ssp[2], "/bio/", "CHELSA_bio1_2071-2100_", lowerc_model, "_", ssp[2], "_V.2.1.tif"))
      cropped_2 <- terra::crop(bio1_fut_2, env_pres[[1]], mask = TRUE) # is this necessary?
      resamp_2 <- terra::resample(cropped_2, env_pres[[1]])

      env_fut <- c(resamp_1, resamp_2, resamp_ndvi)
    }
  }

  if (!future) env_fut <- NULL

  return(list(env_pres = env_pres, env_fut = env_fut))
}

#' Export raw values from RDA model result
#'
#' @param mod RDA model object
#' @param output_path path to all RDA output files
#' @param suffix file suffix
#'
#' @return
#' @export
export_rda_files <- function(mod, output_path, suffix) {
  data.frame(mod$colsum) %>% rownames_to_column(var = "locus") %>% write_csv(file = paste0(output_path, "/RDA_colsum_", suffix, ".csv"), col_names = TRUE)
  data.frame(mod$Ybar) %>% rownames_to_column(var = "INDV") %>% write_csv(file = paste0(output_path, "/RDA_ybar_", suffix, ".csv"), col_names = TRUE)
  data.frame(mod$CCA$v) %>% rownames_to_column(var = "locus") %>% write_csv(file = paste0(output_path, "/RDA_v_", suffix, ".csv"), col_names = TRUE)
  data.frame(mod$CCA$u) %>% rownames_to_column(var = "INDV") %>% write_csv(file = paste0(output_path, "/RDA_u_", suffix, ".csv"), col_names = TRUE)
  data.frame(mod$CCA$wa) %>% rownames_to_column(var = "INDV") %>% write_csv(file = paste0(output_path, "/RDA_wa_", suffix, ".csv"), col_names = TRUE)
  data.frame(mod$CCA$QR$qr) %>% write_csv(file = paste0(output_path, "/RDA_qr_", suffix, ".csv"), col_names = TRUE)
  data.frame(mod$CCA$eig) %>% rownames_to_column(var = "RDA") %>% write_csv(file = paste0(output_path, "/RDA_eig_", suffix, ".csv"), col_names = TRUE)
  data.frame(mod$CCA$biplot) %>% rownames_to_column(var = "var") %>% write_csv(file = paste0(output_path, "/RDA_biplot_", suffix, ".csv"), col_names = TRUE)
  data.frame(mod$CCA$QR$qraux) %>% write_csv(file = paste0(output_path, "/RDA_qraux_", suffix, ".csv"), col_names = TRUE)
  data.frame(mod$CCA$envcentre) %>% tibble::rownames_to_column(var = "axis") %>% write_csv(file = paste0(output_path, "/RDA_envcentre_", suffix, ".csv"), col_names = TRUE)
  data.frame(mod_chi = mod$tot.chi,
            mod_chi_cca = mod$CCA$tot.chi) %>% write_csv(file = paste0(output_path, "/RDA_chi_", suffix, ".csv"), col_names = TRUE)
  # Export scores (scaled and unscaled); labeled "species" corresponds to v, "sites" corresponds to wa, and "constraints" corresponds to u
  scaled_loadings <- vegan::scores(mod, choices = 1:ncol(mod$CCA$v), tidy = TRUE) %>% write_csv(file = paste0(output_path, "/RDA_scaledloadings_", suffix, ".csv"), col_names = TRUE)
  unscaled_loadings <- vegan::scores(mod, choices = 1:ncol(mod$CCA$v), tidy = TRUE, scaling = 0) %>% write_csv(file = paste0(output_path, "/RDA_unscaledloadings_", suffix, ".csv"), col_names = TRUE)
}

#' Import RDA files relevant for adaptive index
#'
#' @param output_path path to all RDA output files
#' @param suffix file suffix
#' @param rds_obj if TRUE, will import only RDS file; defaults to FALSE
#'
#' @return
#' @export
import_rda <- function(output_path, suffix, rds_obj = FALSE) {
  if (!rds_obj) {
    mod <- NULL
    biplot <- readr::read_csv(paste0(output_path, "/RDA_biplot_", suffix, ".csv"), col_names = TRUE) %>% tibble::column_to_rownames(var = "var")
    scaled_loadings <- readr::read_csv(paste0(output_path, "/RDA_scaledloadings_", suffix, ".csv"), col_names = TRUE)
    unscaled_loadings <- readr::read_csv(paste0(output_path, "/RDA_unscaledloadings_", suffix, ".csv"), col_names = TRUE)
    eig <- read_tsv(paste0(output_path, "/RDA_eig_", suffix, ".csv")) %>% column_to_rownames(var = "RDA")
  }
  if (rds_obj) {
    mod <- readRDS(paste0(output_path, "/", suffix, ".RDS"))
    biplot <- data.frame(mod$CCA$biplot)
    eig <- data.frame(mod$CCA$eig)
    # loadings <- as.data.frame(vegan::scores(mod, choices = 1:ncol(mod$CCA$v), display = "species")) %>% rownames_to_column(var = "locus")
    scaled_loadings <- vegan::scores(mod, choices = 1:ncol(mod$CCA$v), tidy = TRUE)
    unscaled_loadings <- vegan::scores(mod, choices = 1:ncol(mod$CCA$v), tidy = TRUE, scaling = 0)
  }
  return(list(mod = mod, biplot = biplot, eig = eig, scaledload = scaled_loadings, unscaledload = unscaled_loadings))
}


# Adaptive index ----------------------------------------------------------

#' Project adaptive component turnover across the landscape
#' Code adapted from Capblancq & Forester (2021) https://doi.org/10.1111/2041-210X.13722
#' GitHub repo available here: https://github.com/Capblancq/RDA-landscape-genomics/blob/main/src/adaptive_index.R
#' TODO convert from raster to terra
#' TODO deal with CRS - currently it's in lat/long
#'
#' @param biplot rownames must match naming of env layers
#' @param K number of PCs to retain if `method = "loadings"` (defaults to 3)
#' @param env_pres environmental layer(s) on which to project
#' @param coords sampling coords (only x and y)
#' @param mod if `method = "predict"` selected, RDA model
#'
#' @return
#' @export
adaptive_index <- function(biplot, K = 3, env_pres, coords, mod = NULL) {
  # Extract values from our environmental rasters
  env <- terra::extract(env_pres, coords)
  # Standardize environmental variables and make into dataframe
  env <- scale(env, center = TRUE, scale = TRUE)
  # Recovering scaling coefficients for extracted env values
  scale_env <- attr(env, 'scaled:scale')
  center_env <- attr(env, 'scaled:center')

  # Formatting environmental rasters for projection
  # var_env_proj_pres <- as.data.frame(raster::rasterToPoints(env_pres[[row.names(biplot)]]))
  var_env_proj_pres <- terra::as.data.frame(env_pres[[row.names(biplot)]], xy = TRUE, na.rm = FALSE)

  # Standardization of the environmental variables based on scaling coefficients
  var_env_proj_RDA <- as.data.frame(scale(var_env_proj_pres[,-c(1,2)],
                                          center_env[row.names(biplot)],
                                          scale_env[row.names(biplot)]))

  # Predicting pixels' genetic component based on RDA axes
  Proj_pres <- list()
  for (i in 1:K) {
    ras_pres <- raster::rasterFromXYZ(data.frame(var_env_proj_pres[,c(1,2)], Z = as.vector(apply(var_env_proj_RDA[,rownames(biplot[i])], 1, function(x) sum(x * biplot[i])))), crs = raster::crs(env_pres))
    # terra::rast(data.frame(var_env_proj_pres[,c(1,2)], type = "xyz")
    names(ras_pres) <- paste0("RDA_pres_", as.character(i))
    Proj_pres[[i]] <- ras_pres
    names(Proj_pres)[i] <- paste0("RDA", as.character(i))
  }

  # Returning projections for current climates for each RDA axis
  return(Proj_pres = Proj_pres)
}


# Plotting AI/offset ------------------------------------------------------

#' Build projected map of RDA results (e.g., RDA adaptive index)
#'
#' @param Proj_data data to plot; either list of rasters, RasterStack, or SpatRaster
#' @param bkg shapefile for plotting (e.g., California state)
#' @param to_mask whether projection needs to be cropped and masked to bkg (defaults to FALSE); only works if Proj_data is raster
#' @param index_name name of index that's being plotted (defaults to "Adaptive index")
#' @param plot_type style for plotting; options are "capblancq" from Capblancq paper, "basic" (default), "rainbow" for rainbow plots with biplot legend, "extracted_vals" for extracted values, or "extracted_rainbow" for extracted values from rainbow map for sampling coordinates
#' @param loadings if `plot_type = "rainbow"` or `plot_type = "extracted_rainbow"`, scaled loadings from RDA model
#' @param coords if `plot_type = "rainbow"` or `plot_type = "extracted_vals"` or `plot_type = "extracted_rainbow"`, sampling coordinates
#' @param biplot_axes if `plot_type = "rainbow"`, which RDA axes to include (defaults to 1 and 2)
#' @param viridis_option option for viridis coloring
#' @param bkg_col if `plot_type = "rainbow"`, background color for blank layer (defaults to "white")
#'
#' @return ggplot2 plot object
#' @export
plot_adaptive <- function(Proj_data, bkg, to_mask = FALSE, index_name = "Adaptive index",
                          plot_type = "basic", loadings, coords, biplot_axes = c(1, 2),
                          viridis_option = "D", bkg_col = "white") {

  if (inherits(Proj_data, "SpatRaster")) Proj_data <- raster::stack(Proj_data)
  n_layers = raster::nlayers(Proj_data)

  # Vectorization of the climatic rasters for ggplot
  RDA_proj <- list()
  for (i in 1:n_layers) {
    RDA_proj[[i]] <- Proj_data[[i]]
  }

  # Turn rasters into points for ggplot
  RDA_proj <- lapply(RDA_proj, function(x) raster::rasterToPoints(x))
  # For each RDA axis, grab value and scale to min and max for that axis
  for (i in 1:length(RDA_proj)) {
    RDA_proj[[i]][,3] <- (RDA_proj[[i]][,3] - min(RDA_proj[[i]][,3]))/(max(RDA_proj[[i]][,3]) - min(RDA_proj[[i]][,3]))
  }

  # Adaptive genetic turnover projected for RDA indexes
  # Bind together points from both RDA axes
  TAB_RDA <- as.data.frame(do.call(rbind, RDA_proj[1:n_layers]))
  colnames(TAB_RDA)[3] <- "value"
  # Add another column to df that specifies which RDA axis values are coming from
  if (n_layers == 1) TAB_RDA$variable <- factor(c(rep("RDA1", nrow(RDA_proj[[1]]))), levels = c("RDA1"))
  if (n_layers == 2) TAB_RDA$variable <- factor(c(rep("RDA1", nrow(RDA_proj[[1]])), rep("RDA2", nrow(RDA_proj[[2]]))), levels = c("RDA1", "RDA2"))
  if (n_layers == 3) TAB_RDA$variable <- factor(c(rep("RDA1", nrow(RDA_proj[[1]])), rep("RDA2", nrow(RDA_proj[[2]])), rep("RDA3", nrow(RDA_proj[[3]]))), levels = c("RDA1", "RDA2", "RDA3"))

  # write_tsv(TAB_RDA, here(output_path, "AI_dat.txt"), col_names = TRUE)

  # Make plot, style is from original Capblancq paper
  if (plot_type == "capblancq") {
    p <- ggplot2::ggplot(data = TAB_RDA) +
      ggplot2::geom_sf(data = bkg, fill = "lightgrey") +
      ggplot2::geom_raster(aes(x = x, y = y, fill = cut(value, breaks = seq(0, 1, length.out = 10), include.lowest = T))) +
      scale_fill_viridis_d(alpha = 0.8, direction = -1, option = "A", labels = c("Negative scores","","","","Intermediate scores","","","","Positive scores")) +
      ggplot2::scale_fill_viridis_d(alpha = 0.8, direction = -1, labels = c("Negative","","","","Intermediate","","","","Positive")) +
      ggplot2::geom_sf(data = bkg, fill = NA, size = 0.1) +
      ggplot2::xlab("Longitude") +
      ggplot2::ylab("Latitude") +
      ggplot2::guides(fill = guide_legend(title = paste0(index_name))) +
      ggplot2::facet_grid(~variable) +
      cowplot::theme_map() +
      ggplot2::theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size = 11))
  }

  if (plot_type == "basic") {
      p <- ggplot2::ggplot(data = TAB_RDA) +
        ggplot2::geom_sf(data = bkg, fill = "lightgrey") +
        ggplot2::geom_raster(aes(x = x, y = y, fill = value)) +
        ggplot2::scale_fill_viridis_c(option = viridis_option, alpha = 0.8, direction = -1) +
        ggplot2::geom_sf(data = bkg, fill = NA, size = 0.1) +
        ggplot2::xlab("Longitude") +
        ggplot2::ylab("Latitude") +
        ggplot2::guides(fill = guide_legend(title = paste0(index_name))) +
        ggplot2::facet_grid(~variable) +
        cowplot::theme_map() +
        ggplot2::theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size = 11))
    }

  if (plot_type == "extracted_vals") {
    rmin <- list()
    rmax <- list()
    rscale <- list()
    rmin[[1]] = minValue(Proj_data[[1]])
    rmax[[1]] = maxValue(Proj_data[[1]])
    rmin[[2]] = minValue(Proj_data[[2]])
    rmax[[2]] = maxValue(Proj_data[[2]])
    rscale[[1]] <- ((Proj_data[[1]] - rmin[[1]]) / (rmax[[1]] - rmin[[1]]))
    rscale[[2]] <- ((Proj_data[[2]] - rmin[[2]]) / (rmax[[2]] - rmin[[2]]))

    if (n_layers == 2) rast <- c(terra::rast(rscale[[1]]), terra::rast(rscale[[2]]))
    if (n_layers == 3) {
      rmin[[3]] = minValue(Proj_data[[3]])
      rmax[[3]] = maxValue(Proj_data[[3]])
      rast <- c(terra::rast(rscale[[1]]), terra::rast(rscale[[2]]), terra::rast(rscale[[3]]))
    }

    ext <- terra::extract(rast, coords, ID = FALSE, xy = TRUE)

    p <- ggplot() +
      geom_sf(data = bkg, color = "lightgrey") +
      geom_point(data = ext %>% pivot_longer(cols = 1:n_layers, names_to = "axis", values_to = "value"),
                 aes(x = x, y = y, color = value), size = 3) +
      scale_color_viridis_c(option = viridis_option, na.value = "transparent", direction = -1, limits = c(0, 1)) +
      facet_grid(~axis) +
      theme_map()
  }

  if (plot_type == "rainbow" | plot_type == "extracted_rainbow") {
    rainbow <- rainbow_map(Proj_data = Proj_data, bkg = bkg, n_layers = n_layers,
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

#' Build rainbow map with vector loadings
#'
#' @param Proj_data
#' @param bkg
#' @param n_layers
#' @param loadings
#' @param biplot_axes
#' @param coords
#' @param bkg_col
#'
#' @returns
#' @export
rainbow_map <- function(Proj_data, bkg, n_layers, loadings, biplot_axes, coords, bkg_col) {
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
    warning("Fewer than three non-zero coefficients provided, adding black/white substitute layers to RGB plot")
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
  
  # Build base map
  p_rainbow <- ggplot() +
    geom_sf(data = bkg, fill = "lightgrey") +
    theme_map()

  if(bkg_col == "white") p_rainbow <- p_rainbow + geom_spatraster_rgb(data = aiRGB, r = 2, g = 1, b = 3)
  if(bkg_col == "black") p_rainbow <- p_rainbow + geom_spatraster_rgb(data = aiRGB, r = 3, g = 2, b = 1)

  p_var <- plot_var_loadings(Proj_data = Proj_data, 
                             loadings = loadings,
                             biplot_axes = biplot_axes, 
                             aiRGB = aiRGB, 
                             coords = coords)
  return(list(map = p_rainbow, vector_load = p_var$var_load_plot, pcacols = p_var$pcacols, rastRGB = aiRGB))
}

#' Plot variable vector loadings as legend for rainbow map
#'
#' @param Proj_data
#' @param loadings
#' @param biplot_axes
#' @param aiRGB
#' @param coords
#'
#' @returns
#' @export
plot_var_loadings <- function(Proj_data, loadings, biplot_axes, aiRGB, coords) {
  TAB_inds <- data.frame(names = rownames(loadings %>% filter(score == "sites")), loadings %>% filter(score == "sites"))
  TAB_var <- data.frame(names = rownames(loadings %>% filter(score == "biplot")), loadings %>% filter(score == "biplot"))

  # Select axes for plotting
  xax <- paste0("RDA", biplot_axes[1])
  yax <- paste0("RDA", biplot_axes[2])
  TAB_inds_sub <- TAB_inds %>% dplyr::select(any_of(c(xax, yax)))
  colnames(TAB_inds_sub) <- c("x", "y")
  TAB_var_sub <- TAB_var %>% column_to_rownames(var = "label") %>% dplyr::select(c(xax, yax))
  colnames(TAB_var_sub) <- c("x", "y")

  # Scale the variable loadings for the arrows
  TAB_var_sub$x <- TAB_var_sub$x * max(TAB_inds_sub$x) / stats::quantile(TAB_var_sub$x)[4]
  TAB_var_sub$y <- TAB_var_sub$y * max(TAB_inds_sub$y) / stats::quantile(TAB_var_sub$y)[4]

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
  plot <- ggplot() +
    geom_hline(yintercept = 0, linewidth = 0.5, col = "gray") +
    geom_vline(xintercept = 0, linewidth = 0.5, col = "gray") +
    geom_point(data = rastvals, aes(x = x, y = y), col = rastpcacols, size = 4, alpha = 0.02) +
    geom_point(data = TAB_inds_sub, aes(x = x, y = y), fill = pcacols, col = "black", pch = 21, size = 3) +
    geom_segment(data = TAB_var_sub, aes(xend = x, yend = y, x = 0, y = 0), color = "black", linewidth = 0.15, linetype = 1, arrow = ggplot2::arrow(length = ggplot2::unit(0.02, "npc"))) +
    ggrepel::geom_text_repel(data = TAB_var_sub, aes(x = x, y = y, label = rownames(TAB_var_sub)), size = 4) +
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

#' Scale a raster stack from 0 to 255
#'
#' @param s RasterStack
#'
#' @noRd
#' @export
stack_to_rgb <- function(s) {
  stack_list <- as.list(s)
  new_stack <- terra::rast(purrr::map(stack_list, raster_to_rgb))
  return(new_stack)
}

#' Scale raster from 0 to 255
#'
#' @param r SpatRast
#'
#' @noRd
#' @export
raster_to_rgb <- function(r) {
  rmax <- terra::minmax(r)["max", ]
  rmin <- terra::minmax(r)["min", ]
  if ((rmax - rmin) == 0) {
    r[] <- 255
  } else {
    r <- (r - rmin) / (rmax - rmin) * 255
  }
  return(r)
}

#' Helper function to create rgb vector
#'
#' @export
#' @noRd
create_rgb_vec <- function(vec) {
  if (any(is.na(vec))) x <- NA else x <- rgb(vec[2], vec[1], vec[3], maxColorValue = 255)
  return(x)
}

# RDA biplot --------------------------------------------------------------

#' Make a biplot of RDA results
#'
#' @param biplot biplot output from RDA model
#' @param loadings (scaled?) loadings from RDA model
#' @param pvals p-values for outlier SNPs
#' @param sig alpha threshold (defaults to 0.01)
#' @param col_w_pvals name of column with p-values to be used
#' @param biplot_type how points are colorized on plot; options are "overall" or "separate"
#' @param biplot_axes if biplot_type = "overall", which RDA axes to plot (defaults to 1 and 2)
#'
#' @return
#' @export
biplot_plot <- function(biplot, loadings, pvals, sig = 0.01, col_w_pvals, biplot_type = "overall", biplot_axes = c(1, 2)) {
  if (biplot_type == "overall") {
    tidy_dat <- tidy_biplot_helper(biplot, biplot_axes, loadings, pvals, sig, col_w_pvals)
    # Biplot of RDA SNPs and scores for variables
    plt_biplot <- biplot_helper(TAB_snps_sub = tidy_dat$TAB_snps_sub,
                                TAB_var_sub = tidy_dat$TAB_var_sub,
                                biplot_type = "overall")
  }
  if (biplot_type == "separate") {
    tidy_dat_12 <- tidy_biplot_helper(biplot, biplot_axes = c(1, 2), loadings, pvals, sig, col_w_pvals)
    tidy_dat_13 <- tidy_biplot_helper(biplot, biplot_axes = c(1, 3), loadings, pvals, sig, col_w_pvals)
    p_rda1 <- biplot_helper(TAB_snps_sub = tidy_dat_12$TAB_snps_sub,
                           TAB_var_sub = tidy_dat_12$TAB_var_sub,
                           biplot_type = "separate",
                           color_by = "x",
                           xax = tidy_dat_12$xax,
                           yax = tidy_dat_12$yax)
    p_rda2 <- biplot_helper(TAB_snps_sub = tidy_dat_12$TAB_snps_sub,
                            TAB_var_sub = tidy_dat_12$TAB_var_sub,
                            biplot_type = "separate",
                            color_by = "y",
                           xax = tidy_dat_12$xax,
                           yax = tidy_dat_12$yax)
    p_rda3 <- biplot_helper(TAB_snps_sub = tidy_dat_13$TAB_snps_sub,
                            TAB_var_sub = tidy_dat_13$TAB_var_sub,
                            biplot_type = "separate",
                            color_by = "y",
                           xax = tidy_dat_13$xax,
                           yax = tidy_dat_13$yax)
    plt_biplot <- list(p_rda1 = p_rda1, p_rda2 = p_rda2, p_rda3 = p_rda3)
  }
  return(plt_biplot)
}

#' Retrieve axis loadings for relevant axes to be plotted
#'
#' @param biplot
#' @param biplot_axes
#' @param loadings
#' @param pvals
#' @param sig
#' @param col_w_pvals
#'
#' @return
#' @export
tidy_biplot_helper <- function(biplot, biplot_axes, loadings, pvals, sig, col_w_pvals) {
  # Select axes for plotting
  xax <- paste0("RDA", biplot_axes[1])
  yax <- paste0("RDA", biplot_axes[2])
  TAB_var <- biplot
  TAB_snps <- ggtidy_helper(loadings, pvals, sig, col_w_pvals)
  TAB_snps_sub <- TAB_snps[, c(xax, yax, "type")]
  colnames(TAB_snps_sub) <- c("x", "y", "type")
  TAB_var_sub <- TAB_var[, c(xax, yax)]
  colnames(TAB_var_sub) <- c("x", "y")

  # Scale the variable loadings for the arrows
  TAB_var_sub$x <- TAB_var_sub$x * max(TAB_snps_sub$x) / stats::quantile(TAB_var_sub$x)[4]
  TAB_var_sub$y <- TAB_var_sub$y * max(TAB_snps_sub$y) / stats::quantile(TAB_var_sub$y)[4]

  return(list(TAB_var_sub = TAB_var_sub, TAB_snps_sub = TAB_snps_sub, xax = xax, yax = yax))
}

ggtidy_helper <- function(loadings, pvals, sig, col_w_pvals) {
  TAB_snps <- loadings %>% dplyr::filter(score == "species") %>% rename(locus = label)
  # TODO fix feeding colname into below
  TAB_df <- full_join(TAB_snps, pvals) %>%
    mutate(pos = 1:nrow(TAB_snps),
           type = case_when(q.values <= sig ~ "Outlier",
                            q.values > sig ~ "Non-outlier")) %>%
    # tidyr::separate_wider_delim(cols = locus,
    #                             delim = "_",
    #                             names = c("chrom", "site", "ref", "alt")
    tidyr::extract(locus, into = c("chrom", "site", "ref", "alt"),
                   regex = "(.*)_([^_]+)_([^_]+)_([^_]+)$")
  return(TAB_df)
}

#' Helper function to build RDA biplot
#'
#' @param TAB_snps_sub df with SNPs and loadings
#' @param TAB_var_sub df with variable loadings
#' @param biplot_type "overall" or "separate" if you want SNPs colorized by axis
#' @param color_by if biplot_type = "separate", which axis to color SNPs by ("x" or "y")
#' @param xax label for x axis
#' @param yax label for y axis
#'
#' @return
#' @export
biplot_helper <- function(TAB_snps_sub, TAB_var_sub, xax, yax, biplot_type, color_by) {
  if (biplot_type == "overall") {
    plt_biplot <-
    ggplot2::ggplot() +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = gray(.80), linewidth = 0.6) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = gray(.80), linewidth = 0.6) +
      ggplot2::geom_point(data = TAB_snps_sub, ggplot2::aes(x = x, y = y, colour = type), size = 1.4) +
      ggplot2::scale_color_manual(values = c(rgb(0.7, 0.7, 0.7, 0.1), "#F9A242FF")) +
      ggplot2::geom_segment(data = TAB_var_sub, ggplot2::aes(xend = x, yend = y, x = 0, y = 0), colour = "black", linewidth = 0.5, linetype = 1, arrow = ggplot2::arrow(length = ggplot2::unit(0.02, "npc"))) +
      ggrepel::geom_text_repel(data = TAB_var_sub, ggplot2::aes(x = x, y = y, label = row.names(TAB_var_sub)), size = 4) +
      ggplot2::xlab(xax) +
      ggplot2::ylab(yax) +
      ggplot2::guides(color = ggplot2::guide_legend(title = "SNP type")) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(panel.background = ggplot2::element_blank(),
                     legend.background = ggplot2::element_blank(),
                     panel.grid = ggplot2::element_blank(),
                     plot.background = ggplot2::element_blank(),
                     legend.text = ggplot2::element_text(size = ggplot2::rel(.8)),
                     strip.text = ggplot2::element_text(size = 11))
  }
  if (biplot_type == "separate") {
    if (color_by == "x") {
      plt_biplot <-
      ggplot2::ggplot() +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = gray(.80), linewidth = 0.6) +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = gray(.80), linewidth = 0.6) +
        ggplot2::geom_point(data = TAB_snps_sub %>% filter(type == "Non-outlier"), ggplot2::aes(x = x, y = y), size = 1.4, alpha = 0.1, color = "grey") +
        ggplot2::geom_point(data = TAB_snps_sub %>% filter(type == "Outlier"), ggplot2::aes(x = x, y = y, colour = x), size = 1.4) +
        scale_color_viridis_c(alpha = 0.8, direction = -1, name = xax) +
        ggplot2::geom_segment(data = TAB_var_sub, ggplot2::aes(xend = x, yend = y, x = 0, y = 0), colour = "black", linewidth = 0.5, linetype = 1, arrow = ggplot2::arrow(length = ggplot2::unit(0.02, "npc"))) +
        ggrepel::geom_text_repel(data = TAB_var_sub, ggplot2::aes(x = x, y = y, label = row.names(TAB_var_sub)), size = 4) +
        ggplot2::xlab(xax) +
        ggplot2::ylab(yax) +
        ggplot2::theme_bw(base_size = 11) +
        ggplot2::theme(panel.background = ggplot2::element_blank(),
                       # legend.background = ggplot2::element_blank(),
                       panel.grid = ggplot2::element_blank(),
                       plot.background = ggplot2::element_blank(),
                       # legend.text = ggplot2::element_text(size = ggplot2::rel(.8)),
                       legend.position = "none",
                       strip.text = ggplot2::element_text(size = 11))
    }
    if (color_by == "y") {
      plt_biplot <-
      ggplot2::ggplot() +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = gray(.80), linewidth = 0.6) +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = gray(.80), linewidth = 0.6) +
        ggplot2::geom_point(data = TAB_snps_sub %>% filter(type == "Non-outlier"), ggplot2::aes(x = x, y = y), size = 1.4, alpha = 0.1, color = "grey") +
        ggplot2::geom_point(data = TAB_snps_sub %>% filter(type == "Outlier"), ggplot2::aes(x = x, y = y, colour = y), size = 1.4) +
        scale_color_viridis_c(alpha = 0.8, direction = -1, name = yax) +
        ggplot2::geom_segment(data = TAB_var_sub, ggplot2::aes(xend = x, yend = y, x = 0, y = 0), colour = "black", linewidth = 0.5, linetype = 1, arrow = ggplot2::arrow(length = ggplot2::unit(0.02, "npc"))) +
        ggrepel::geom_text_repel(data = TAB_var_sub, ggplot2::aes(x = x, y = y, label = row.names(TAB_var_sub)), size = 4) +
        ggplot2::xlab(xax) +
        ggplot2::ylab(yax) +
        ggplot2::theme_bw(base_size = 11) +
        ggplot2::theme(panel.background = ggplot2::element_blank(),
                       # legend.background = ggplot2::element_blank(),
                       panel.grid = ggplot2::element_blank(),
                       plot.background = ggplot2::element_blank(),
                       legend.position = "none",
                       # legend.text = ggplot2::element_text(size = ggplot2::rel(.8)),
                       strip.text = ggplot2::element_text(size = 11))
    }

  }
  return(plt_biplot)
}
