krig_raster_auto <- function(r, weight_r = NULL,
                             candidate_models = c("Sph", "Exp", "Gau", "Mat"),
                             max_range_frac = 0.5, nmax = 30, verbose = TRUE) {
  # --- Step 1: Convert raster to points ---
  r_pts <- terra::as.data.frame(r, xy = TRUE, na.rm = TRUE)
  names(r_pts)[3] <- "value"
  
  # Add weights if weight_r is provided
  if (!is.null(weight_r)) {
    weights <- terra::extract(weight_r, r_pts[, c("x", "y")])[,2]
    r_pts$weight <- weights
    r_pts$weight[is.na(r_pts$weight) | r_pts$weight <= 0] <- 1  # Avoid NA/zero weights
    # Scale weights: higher sample counts → higher weight
    r_pts$weight <- r_pts$weight / max(r_pts$weight, na.rm = TRUE)
  } else {
    r_pts$weight <- NULL  # No weights
  }
  
  # Convert to sf object
  r_sf <- sf::st_as_sf(r_pts, coords = c("x", "y"), crs = terra::crs(r))
  
  # --- Step 2: Fit variogram (unweighted) ---
  if (verbose) cat("Fitting empirical variogram (unweighted)...\n")
  formula_str <- value ~ 1
  v <- gstat::variogram(formula_str, data = r_sf)
  
  if (verbose) cat("Fitting theoretical variogram models...\n")
  # Approximate starting values
  psill_start <- stats::var(r_pts$value, na.rm = TRUE) * 0.8
  nugget_start <- stats::var(r_pts$value, na.rm = TRUE) * 0.2
  range_start <- max(v$dist, na.rm = TRUE) * max_range_frac
  
  # Fit multiple models and pick the best
  fitted_models <- lapply(candidate_models, function(mod) {
    tryCatch({
      start_model <- gstat::vgm(psill = psill_start, model = mod, range = range_start, nugget = nugget_start)
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
  
  # --- Step 3: Kriging (apply weights here) ---
  if (verbose) cat("Performing Kriging prediction...\n")
  krig_model <- gstat::gstat(formula = value ~ 1, locations = r_sf,
                             model = best_fit, weights = r_pts$weight, nmax = nmax)
  
  # Create prediction grid
  grid <- terra::as.data.frame(r, xy = TRUE, na.rm = FALSE)
  grid_sf <- sf::st_as_sf(grid, coords = c("x", "y"), crs = terra::crs(r))
  krig_pred <- predict(krig_model, newdata = grid_sf)
  
  # --- Step 4: Convert predictions to rasters ---
  krig_df <- as.data.frame(krig_pred)
  r_pred <- terra::rast(terra::ext(r), resolution = terra::res(r), crs = terra::crs(r))
  terra::values(r_pred) <- krig_df$var1.pred
  
  r_var <- terra::rast(terra::ext(r), resolution = terra::res(r), crs = terra::crs(r))
  terra::values(r_var) <- krig_df$var1.var
  
  return(list(prediction = r_pred, variance = r_var, variogram = v, model = best_fit))
}


r <- rast(here("analysis/genetic_diversity/outputs/", "wingen.tif"))[[1]] 
weight_r <-  rast(here("analysis/genetic_diversity/outputs/", "wingen.tif"))[[2]]  # Replace with your weight raster path
results <- krig_raster_auto(
  r = r,
  weight_r = weight_r,
  verbose = TRUE
)

kriged_r <- results$prediction

ca_parts <- st_cast(ca_proj, "POLYGON")
ca_parts$area <- st_area(ca_parts)
ca_mainland <- ca_parts[which.max(ca_parts$area), ]
weighted_krig <- mask(mask(kriged_r, range_map), ca_mainland)

terra::plot(results$prediction, main = "Kriging Prediction (Higher Samples = Higher Weight)")
terra::plot(results$variance, main = "Kriging Variance")










# OLD
# --- Step 1: Load raster and convert to points ---
# Replace with your raster path
r <- rast(here("analysis/genetic_diversity/outputs/", "wingen.tif"))[[1]]  # Or use a built-in raster

# Convert raster to data frame (x, y, value)
r_pts <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
names(r_pts)[3] <- "value"

# Convert to spatial object for gstat
coordinates(r_pts) <- ~x+y

# Set CRS if missing
crs(r_pts) <- crs(r)

# --- Step 2: Create interpolation grid ---
# Create a regular grid covering the raster extent
grid <- as.data.frame(lyr, xy = TRUE, na.rm = FALSE)  # Use the raster extent and resolution
names(grid) <- c("x", "y", "value")
coordinates(grid) <- ~x+y
crs(grid) <- crs(r)

# --- Step 3: Variogram and model fitting ---
# Empirical variogram
library(gstat)
v <- variogram(value ~ 1, data = r_pts)

# Fit a variogram model (exponential as default)
v_fit <- fit.variogram(v, model = vgm(psill = var(r_pts$value, na.rm = TRUE), model = "Exp", range = 3.5e5, nugget = 1e-08))

# Plot variogram and fitted model
pdf(here("analysis", "genetic_diversity", "plots", "wingen_variogram.pdf"))
plot(v, v_fit, main = "Empirical Variogram & Fitted Model")
dev.off()

# --- Step 4: Kriging prediction ---
# Set up Kriging model (ordinary Kriging)
krig_model <- gstat(formula = value ~ 1, locations = r_pts, model = v_fit, nmax = 30)  # nmax limits neighbors for speed

# Predict on grid
krig_pred <- predict(krig_model, newdata = grid)

# --- Step 5: Convert prediction to raster ---
# Back to data frame
krig_df <- as.data.frame(krig_pred)

# Create prediction raster
r_krig <- rast(ext(r), resolution = res(r), crs = crs(r))
values(r_krig) <- krig_df$var1.pred

# Create variance raster (uncertainty)
r_var <- rast(ext(r), resolution = res(r), crs = crs(r))
values(r_var) <- krig_df$var1.var

unweighted_krig <- mask(mask(r_krig, range_map), ca_proj)  # Mask to remove low variance areas

# --- Step 6: Plot results ---
pdf(here("analysis", "genetic_diversity", "plots", "wingen_kriging.pdf"))

plot(r_mask, main = "Kriging Prediction with Mask")
plot(idw_rast, main = "IDW Interpolation")
plot(r_krig, main = "Kriging Prediction")
plot(r_var, main = "Kriging Variance (Uncertainty)")
dev.off()
