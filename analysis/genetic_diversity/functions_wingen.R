#' GAM-Based Interpolation of Windowed Diversity
#'
#' Perform generalized additive model based interpolation of the raster(s) produced by \link[wingen]{window_gd} using \link[mgcv]{gam}.
#'
#' @param r SpatRaster produced by \link[wingen]{window_gd}
#' @param grd object to create grid for interpolation; can be a SpatRaster or RasterLayer. If undefined, will use \code{r} to create a grid.
#' @param index integer index of the layer in the raster stack to interpolate (defaults to 1; i.e., the first layer)
#' @param coords if provided, interpolation will occur based only on values at these coordinates. Can be provided as an sf points, a two-column matrix, or a data.frame representing x and y coordinates
#' @param agg_grd factor to use for aggregation of `grd`, if provided (this will decrease the resolution of the final interpolated raster; defaults to NULL)
#' @param disagg_grd factor to use for disaggregation of `grd`, if provided (this will increase the resolution of the final interpolated raster; defaults to NULL)
#' @param agg_r factor to use for aggregation of `r`, if provided (this will decrease the number of points used in the interpolation model; defaults to NULL)
#' @param disagg_r factor to use for disaggregation, of `r` if provided (this will increase the number of points used in the interpolation model; defaults to NULL)
#' @param lower_bound if TRUE (default), converts all values in the interpolated raster less than the minimum value of the input raster, to that minimum.
#' @param upper_bound if TRUE (default), converts all values in the interpolated raster greater than the maximum value of the input raster, to that maximum.
#' @param resample whether to resample `grd` or `r`. Set to `"r"` to resample `r` to `grd`. Set to `"grd"` to resample `grd` to `r` (defaults to FALSE for no resampling)
#' @param resample_first if aggregation or disaggregation is used in addition to resampling, specifies whether to resample before (resample_first = TRUE) or after (resample_first = FALSE) aggregation/disaggregation (defaults to TRUE)
#' @param bs Character; the basis type for the spatial smoother (passed to \code{mgcv::s()}, default = "tp").
#' @param k Integer; the basis dimension (maximum wiggliness + 1) for the smoother. Defaults to the closest integar value of the square root of the number of values in `r` (or `coords`, if provided). If the number of values is less than 10, defaults to k = 3.
#' @param method Character; smoothing-parameter estimation method for \code{mgcv::gam()} (default = "REML").
#' @param gamma Numeric; inflation factor for the smoothing selection criterion (default = 1).
#' @param select Logical; whether to apply shrinkage to smooth terms (default = TRUE).
#'
#' @return A SpatRaster of the interpolated surface, with the same CRS as `r`.
#'
#' @examples
#' load_mini_ex()
#' wpi <- window_gd(mini_vcf, mini_coords, mini_lyr, L = 10, rarify = TRUE)
#' gampi <- gam_gd(wpi, mini_lyr)
#' plot_gd(gampi, main = "GAM Pi")
#'
#' @export
gam_gd <- function(r, grd = NULL, index = 1, coords = NULL,
                   agg_grd = NULL, disagg_grd = NULL, agg_r = NULL, disagg_r = NULL,
                   lower_bound = TRUE, upper_bound = TRUE,
                   resample = FALSE, resample_first = TRUE,
                   bs = "tp", k = NULL, method = "REML", gamma = 1, select = TRUE) {

  # Ensure terra objects
  if (!inherits(r, "SpatRaster")) r <- terra::rast(r)
  if (!is.null(grd) && inherits(grd, "RasterLayer")) grd <- terra::rast(grd)

  # Subset layer
  if (terra::nlyr(r) > 1) r <- r[[index]]

  # Transform rasters (agg/disagg/resample)
  stk <- raster_transform(
    r               = r,
    grd             = grd,
    agg_grd         = agg_grd,
    disagg_grd      = disagg_grd,
    agg_r           = agg_r,
    disagg_r        = disagg_r,
    resample        = resample,
    resample_first  = resample_first
  )
  r_t <- stk[[names(r)]]
  grd_t <- stk[["grd"]]

  # Build point data.frame
  if (!is.null(coords)) {
    pts_df <- terra::as.data.frame(
      terra::extract(r_t, coords, ID = FALSE, xy = TRUE),
      xy = TRUE
    )
    names(pts_df)[3] <- "layer"
  } else {
    pts_df <- terra::as.data.frame(r_t, xy = TRUE, na.rm = TRUE)
    names(pts_df)[3] <- "layer"
  }

  # Set K if not provided
  if (is.null(k)) {
    n_obs <- nrow(pts_df)
    if (n_obs < 10) {
      k <- 3  # Minimum basis dimension
    } else {
      k <- floor(sqrt(n_obs))  # Default based on number of observations
    }
    message("Setting k = ", k)
  }

  # Fit GAM
  gam_mod <- mgcv::gam(
    formula = layer ~ s(x, y, bs = bs, k = k),
    data    = pts_df,
    method  = method,
    gamma   = gamma,
    select  = select
  )
  
  # Print gam check results
  message("GAM model check results:")
  mgcv::gam.check(gam_mod)

  # Prediction grid
  grd_df <- terra::as.data.frame(grd_t, xy = TRUE, na.rm = FALSE)

  # Predict values
  pred_vals <- stats::predict(gam_mod, newdata = grd_df)
  grd_df$layer <- pred_vals

  # Apply bounds
  if (is.numeric(lower_bound)) {
    grd_df$layer[grd_df$layer < lower_bound] <- lower_bound
  }
  if (is.numeric(upper_bound)) {
    grd_df$layer[grd_df$layer > upper_bound] <- upper_bound
  }
  if (isTRUE(lower_bound)) {
    mn <- min(pts_df$layer, na.rm = TRUE)
    grd_df$layer[grd_df$layer < mn] <- mn
  }
  if (isTRUE(upper_bound)) {
    mx <- max(pts_df$layer, na.rm = TRUE)
    grd_df$layer[grd_df$layer > mx] <- mx
  }

  # Reconstruct SpatRaster
  out <- terra::rast(grd_df[, c("x", "y", "layer")],
    type = "xyz",
    crs  = terra::crs(r_t)
  )
  names(out) <- names(r_t)

  return(out)
}


#' Transform raster
#'
#' @inheritParams krig_gd
#'
#' @noRd
raster_transform <- function(r, grd, resample = FALSE, agg_grd = NULL, disagg_grd = NULL, agg_r = NULL, disagg_r = NULL, resample_first = TRUE) {
  if (terra::nlyr(r) > 1) stop(">1 layer provided for r")
  if (terra::nlyr(grd) > 1) stop(">1 layer provided for grd")

  if (resample_first) {
    if (resample == "r") r <- terra::resample(r, grd)
    if (resample == "grd") grd <- terra::resample(grd, r)
  }

  if (!is.null(agg_grd) & !is.null(disagg_grd)) stop("Both agg_grd and disagg_grd provided, when only one should be provided")
  if (!is.null(agg_grd)) grd <- terra::aggregate(grd, agg_grd)
  if (!is.null(disagg_grd)) grd <- terra::disagg(grd, disagg_grd)

  if (!is.null(agg_r) & !is.null(disagg_r)) stop("Both agg_r and disagg_r provided, when only one should be provided")
  if (!is.null(agg_r)) r <- terra::aggregate(r, agg_r)
  if (!is.null(disagg_r)) r <- terra::disagg(r, disagg_r)

  if (!resample_first) {
    if (resample == "r") r <- terra::resample(r, grd)
    if (resample == "grd") grd <- terra::resample(grd, r)
  }

  s <- list(r, grd)
  names(s) <- c(names(r), "grd")

  return(s)
}

idw_gd <- function(r){
  # Convert raster cells to points (to get x, y, value)
  r = mask_gd(window_Ho[[1]], window_Ho[[2]], minval = 2)

  # Turn raster into points
  r_pts <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
  names(r_pts)[3] <- "value"

  # Create gstat object for IDW
  idw_gstat <- gstat::gstat(formula = value ~ 1, locations = ~x+y, data = r_pts, nmax = 12, set = list(idp = 2))

  # Create prediction grid (same extent/resolution as original raster)
  grid <- as.data.frame(lyr, xy = TRUE, na.rm = FALSE)
  names(grid) <- c("x", "y", "value")

  # Predict using IDW
  idw_pred <- predict(idw_gstat, newdata = grid)

  # Convert result to raster
  idw_rast <- rast(ext(r), resolution = res(r), crs = crs(r))
  values(idw_rast) <- idw_pred$var1.pred

  return(idw_rast)
}

# gstat vignette (Pebesma, 2004) uses nmax = 20.
# Diggle & Ribeiro (2007), Model-based Geostatistics, discuss ~30 neighbors as typical.
# ArcGIS kriging defaults to 12–32 neighbors.
krig_gd2 <- function(r, weight_r = NULL,
                     candidate_models = c("Sph", "Exp", "Gau", "Mat"),
                     max_range_frac = 0.5, nmax = 30, verbose = TRUE,
                     psill_start = NULL, nugget_start = NULL, range_start = NULL) {

  # Convert raster to points
  r_pts <- terra::as.data.frame(r, xy = TRUE, na.rm = TRUE)
  names(r_pts)[3] <- "value"

  # Convert to sf object
  r_sf <- sf::st_as_sf(r_pts, coords = c("x", "y"), crs = terra::crs(r))

  # Handle weights (convert sample counts to variance estimates)
  # Approximate location-specific measurement variance as σ²_global / n:
  # - σ²_global = variance of all raster values (assumes homoskedasticity)
  # - n = sample count in each cell
  # - Var(mean) = σ² / n because averaging reduces variance by factor n
  # Gstat expects weights = 1/variance to reflect confidence in observations.
  # More samples → lower variance → higher weight.
  if (!is.null(weight_r)) {
    sample_counts <- terra::extract(weight_r, r_pts[, c("x", "y")])[,2]
    
    # Confirm no sample count values are less than 1 or NA
    stopifnot(all(sample_counts >= 1))
    stopifnot(all(!is.na(sample_counts))) 

    # Estimate global variance
    sigma2_global <- stats::var(r_pts$value, na.rm = TRUE)

    # Estimate location-specific variance
    r_pts$variance <- sigma2_global / sample_counts

    # Convert variance to weights
    r_pts$weight <- 1 / r_pts$variance  
  } else {
    r_pts$weight <- NULL  # No location-specific variance information
  }

  # --- Step 2: Fit variogram (unweighted) ---
  formula_str <- value ~ 1
  v <- gstat::variogram(formula_str, data = r_sf)

  # Approximate starting values
  if (is.null(psill_start)) psill_start <- stats::var(r_pts$value, na.rm = TRUE) * 0.8
  if (is.null(nugget_start)) nugget_start <- stats::var(r_pts$value, na.rm = TRUE) * 0.2
  if (is.null(range_start)) range_start <- max(v$dist, na.rm = TRUE) * max_range_frac

  # Fit multiple models and pick the best
  fitted_models <- lapply(candidate_models, function(mod) {
    tryCatch({
      start_model <- gstat::vgm(psill = psill_start, model = mod,
                                range = range_start, nugget = nugget_start)
      fit <- gstat::fit.variogram(v, model = start_model, fit.method = 6)
      attr(fit, "model_name") <- mod
      return(fit)
    }, error = function(e) NULL)
  })
  
  fitted_models <- Filter(Negate(is.null), fitted_models)
  if (length(fitted_models) == 0) stop("Variogram fitting failed for all models.")
  
  sse <- sapply(fitted_models, function(f) attr(f, "SSErr"))
  best_fit <- fitted_models[[which.min(sse)]]
  best_model <- attr(best_fit, "model_name")
  
  if (verbose) {
    cat("Best model:", best_model, "\n")
    cat("SSErr:", min(sse), "\n")
  }

  # --- Step 3: Kriging with weights ---
  krig_model <- gstat::gstat(formula = value ~ 1, locations = r_sf,
                             model = best_fit, weights = r_pts$weight, nmax = nmax)

  # Create prediction grid
  grid <- terra::as.data.frame(r, xy = TRUE, na.rm = FALSE)
  grid_sf <- sf::st_as_sf(grid, coords = c("x", "y"), crs = terra::crs(r))
  krig_pred <- predict(krig_model, newdata = grid_sf)

  # Convert predictions to rasters
  krig_df <- as.data.frame(krig_pred)
  r_pred <- terra::rast(terra::ext(r), resolution = terra::res(r), crs = terra::crs(r))
  terra::values(r_pred) <- krig_df$var1.pred

  r_var <- terra::rast(terra::ext(r), resolution = terra::res(r), crs = terra::crs(r))
  terra::values(r_var) <- krig_df$var1.var

  return(list(prediction = r_pred, variance = r_var, variogram = v, model = best_fit))
}



cokrig_gd2 <- function(r, covariates = NULL, weight_r = NULL,
                       candidate_models = c("Sph", "Exp", "Gau", "Mat"),
                       max_range_frac = 0.5, nmax = 30, verbose = TRUE,
                       psill_start = NULL, nugget_start = NULL, range_start = NULL) {

  # --- Step 1: Prepare data ---
  # Convert primary raster to points
  r_pts <- terra::as.data.frame(r, xy = TRUE, na.rm = TRUE)
  names(r_pts)[3] <- "value"

  # Add covariates (handles SpatRaster stack directly)
  if (!is.null(covariates)) {
    cov_vals <- terra::extract(covariates, r_pts[, c("x", "y")])
    cov_vals <- as_tibble(cov_vals[, -1])  # Drop ID column
    names(cov_vals) <- names(covariates)  # Preserve layer names
    r_pts <- bind_cols(r_pts, cov_vals)
  }

  # Convert to sf object
  r_sf <- st_as_sf(r_pts, coords = c("x", "y"), crs = crs(r))

  # --- Step 2: Handle weights ---
  if (!is.null(weight_r)) {
    sample_counts <- terra::extract(weight_r, r_pts[, c("x", "y")])[,2]
    stopifnot(all(sample_counts >= 1), all(!is.na(sample_counts)))
    sigma2_global <- var(r_pts$value, na.rm = TRUE)
    r_pts$variance <- sigma2_global / sample_counts
    r_pts$weight <- 1 / r_pts$variance
  } else {
    r_pts$weight <- NULL
  }

  # --- Step 3: Fit variogram ---
  # Create formula: value ~ covariate1 + covariate2 + ...
  if (!is.null(covariates)) {
    cov_names <- names(covariates)
    formula_str <- as.formula(paste("value ~", paste(cov_names, collapse = " + ")))
  } else {
    formula_str <- value ~ 1
  }

  # Fit variogram
  r_sf <- drop_na(r_sf)
  v <- variogram(formula_str, data = r_sf)
  if (is.null(psill_start)) psill_start <- var(r_pts$value, na.rm = TRUE) * 0.8
  if (is.null(nugget_start)) nugget_start <- var(r_pts$value, na.rm = TRUE) * 0.2
  if (is.null(range_start)) range_start <- max(v$dist, na.rm = TRUE) * max_range_frac

  fitted_models <- lapply(candidate_models, function(mod) {
    tryCatch({
      start_model <- vgm(psill_start, model = mod,
                         range = range_start, nugget = nugget_start)
      fit <- fit.variogram(v, model = start_model, fit.method = 6)
      attr(fit, "model_name") <- mod
      fit
    }, error = function(e) NULL)
  })
  fitted_models <- Filter(Negate(is.null), fitted_models)
  if (length(fitted_models) == 0) stop("Variogram fitting failed.")
  sse <- sapply(fitted_models, function(f) attr(f, "SSErr"))
  best_fit <- fitted_models[[which.min(sse)]]
  best_model <- attr(best_fit, "model_name")

  if (verbose) {
    cat("Best variogram model:", best_model, "SSE =", min(sse), "\n")
  }

  # --- Step 4: Create gstat model ---
  krig_model <- gstat(formula = formula_str, locations = r_sf,
                      model = best_fit, weights = r_pts$weight, nmax = nmax)

  # --- Step 5: Prepare prediction grid ---
  grid <- terra::as.data.frame(r, xy = TRUE, na.rm = FALSE)[, c("x", "y")]
  if (!is.null(covariates)) {
    grid_covs <- terra::extract(covariates, grid[, c("x", "y")])
    grid_covs <- as_tibble(grid_covs[, -1])
    names(grid_covs) <- names(covariates)
    grid <- bind_cols(grid, grid_covs)
  }
  grid_sf <- st_as_sf(grid, coords = c("x", "y"), crs = crs(r))

  # --- Step 6: Predict ---
  pred <- predict(krig_model, newdata = grid_sf)
  
  # --- Step 7: Convert predictions to rasters ---
  grid_sf$prediction <- pred$var1.pred
  grid_sf$variance <- pred$var1.var
  r_pred <- rasterize(grid_sf, r, field = "prediction")
  r_var <- rasterize(grid_sf, r, field = "variance")

  names(r_pred) <- "prediction"
  names(r_var) <- "variance"

  return(list(prediction = r_pred, variance = r_var, variogram = v, model = best_fit))
}
