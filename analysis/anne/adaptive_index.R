#' Project adaptive component turnover across the landscape
#' Code adapted from Capblancq & Forester (2021) https://doi.org/10.1111/2041-210X.13722
#' GitHub repo available here: https://github.com/Capblancq/RDA-landscape-genomics/blob/main/src/adaptive_index.R
#' TODO convert from raster to terra
#' TODO deal with CRS - currently it's in lat/long
#'
#' @param biplot rownames must match naming of env layers
#' @param K number of PCs to retain if `method = "loadings"`
#' @param env_pres environmental layer(s) on which to project
#' @param coords sampling coords (only x and y)
#' @param range whether to mask to species range (defaults to NULL)
#' @param method how to predict values; options are "loadings" (use RDA results; default) or "predict" (use RDA predict)
#'
#' @return
#' @export
adaptive_index <- function(biplot, K, env_pres, coords, range = NULL, method = "loadings"){
  # Extract values from our environmental rasters
  env <- raster::extract(env_pres, coords)
  # Standardize environmental variables and make into dataframe
  env <- scale(env, center = TRUE, scale = TRUE)
  # Recovering scaling coefficients for extracted env values
  scale_env <- attr(env, 'scaled:scale')
  center_env <- attr(env, 'scaled:center')

  # Formatting environmental rasters for projection
  var_env_proj_pres <- as.data.frame(raster::rasterToPoints(env_pres[[row.names(biplot)]]))

  # Standardization of the environmental variables
  var_env_proj_RDA <- as.data.frame(scale(var_env_proj_pres[,-c(1,2)],
                                          center_env[row.names(biplot)],
                                          scale_env[row.names(biplot)]))

  # Predicting pixels genetic component based on RDA axes
  Proj_pres <- list()
  if (method == "loadings") {
    for (i in 1:K) {
      ras_pres <- raster::rasterFromXYZ(data.frame(var_env_proj_pres[,c(1,2)], Z = as.vector(apply(var_env_proj_RDA[,rownames(biplot[i])], 1, function(x) sum(x * biplot[i])))), crs = crs(env_pres))
      names(ras_pres) <- paste0("RDA_pres_", as.character(i))
      Proj_pres[[i]] <- ras_pres
      names(Proj_pres)[i] <- paste0("RDA", as.character(i))
    }
  }

  # Prediction with RDA model and linear combinations
  if (method == "predict") {
    pred <- predict(RDA, var_env_proj_RDA[,rownames(biplot[i])], type = "lc")
    for (i in 1:K) {
      ras_pres <- raster::rasterFromXYZ(data.frame(var_env_proj_pres[,c(1,2)], Z = as.vector(pred[,i])), crs = crs(env_pres))
      names(ras_pres) <- paste0("RDA_pres_", as.character(i))
      Proj_pres[[i]] <- ras_pres
      names(Proj_pres)[i] <- paste0("RDA", as.character(i))
    }
  }

  # Mask with the range if supplied
  if(!is.null(range)){
    Proj_pres <- lapply(Proj_pres, function(x) terra::mask(x, range))
  }

  # Returning projections for current climates for each RDA axis
  return(Proj_pres = Proj_pres)
}

#' Build projected map of RDA results (e.g., RDA adaptive index)
#'
#' @param Proj_data data to plot
#' @param title scaffold name, for labeling plot
#' @param bkg shapefile for plotting (e.g., California state)
#' @param to_mask whether projection needs to be cropped and masked to bkg (defaults to FALSE)
#'
#' @return ggplot2 plot object
#' @export
plot_adaptiveindex <- function(Proj_data, scaf_name, bkg, to_mask = FALSE) {
  if (to_mask) {
    Proj_data$RDA1 <- crop(terra::rast(Proj_data$RDA1), bkg, mask = TRUE)
    Proj_data$RDA2 <- crop(terra::rast(Proj_data$RDA2), bkg, mask = TRUE)
    Proj_data$RDA3 <- crop(terra::rast(Proj_data$RDA3), bkg, mask = TRUE)
    # Convert back to raster objects
    Proj_data$RDA1 <- raster::raster(Proj_data$RDA1)
    Proj_data$RDA2 <- raster::raster(Proj_data$RDA2)
    Proj_data$RDA3 <- raster::raster(Proj_data$RDA3)
  }
  # Vectorization of the climatic rasters for ggplot
  RDA_proj <- list(Proj_data$RDA1, Proj_data$RDA2, Proj_data$RDA3)
  # Turn rasters into points for ggplot
  RDA_proj <- lapply(RDA_proj, function(x) raster::rasterToPoints(x))
  for (i in 1:length(RDA_proj)) {
    RDA_proj[[i]][,3] <- (RDA_proj[[i]][,3] - min(RDA_proj[[i]][,3]))/(max(RDA_proj[[i]][,3]) - min(RDA_proj[[i]][,3]))
  }

  # Adaptive genetic turnover projected for RDA1 and RDA2 indexes
  # Bind together points from both RDA axes
  TAB_RDA <- as.data.frame(do.call(rbind, RDA_proj[1:3]))
  colnames(TAB_RDA)[3] <- "value"
  # Add another column to df that specifies whether values are coming from RDA axes 1 or 2
  TAB_RDA$variable <- factor(c(rep("RDA1", nrow(RDA_proj[[1]])), rep("RDA2", nrow(RDA_proj[[2]])), rep("RDA3", nrow(RDA_proj[[3]]))),
                             levels = c("RDA1", "RDA2", "RDA3"))

  # Make plot
  ggplot(data = TAB_RDA) +
    # geom_sf(data = admin, fill=gray(.9), size = 0) +
    geom_sf(data = bkg, fill = "lightgrey") +
    geom_raster(aes(x = x, y = y, fill = cut(value, breaks = seq(0, 1, length.out = 10), include.lowest = T))) +
    # scale_fill_viridis_d(alpha = 0.8, direction = -1, option = "A", labels = c("Negative scores","","","","Intermediate scores","","","","Positive scores")) +
    scale_fill_viridis_d(alpha = 0.8, direction = -1, labels = c("Negative","","","","Intermediate","","","","Positive")) +
    # scale_fill_viridis_d(alpha = 0.8, direction = -1) +
    geom_sf(data = bkg, fill = NA, size = 0.1) +
    xlab("Longitude") +
    ylab("Latitude") +
    guides(fill = guide_legend(title = "Adaptive index")) +
    facet_grid(~variable) +
    # theme_bw(base_size = 11) +
    theme_map() +
    theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size = 11)) +
    ggtitle(paste0(scaf_name))
}

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
                           color_by = "x")
    p_rda2 <- biplot_helper(TAB_snps_sub = tidy_dat_12$TAB_snps_sub,
                            TAB_var_sub = tidy_dat_12$TAB_var_sub,
                            biplot_type = "separate",
                            color_by = "y")
    p_rda3 <- biplot_helper(TAB_snps_sub = tidy_dat_13$TAB_snps_sub,
                            TAB_var_sub = tidy_dat_13$TAB_var_sub,
                            biplot_type = "separate",
                            color_by = "y")
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

  return(list(TAB_var_sub = TAB_var_sub, TAB_snps_sub = TAB_snps_sub))
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
#'
#' @return
#' @export
biplot_helper <- function(TAB_snps_sub, TAB_var_sub, biplot_type, color_by) {
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
                     strip.text = ggplot2::element_text(size = 11)) +
      coord_equal()
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
                       strip.text = ggplot2::element_text(size = 11)) +
        coord_equal()
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
                       strip.text = ggplot2::element_text(size = 11)) +
        coord_equal()
    }

  }
  return(plt_biplot)
}
