
library(terra)
library(here)
library(tidyverse)
library(raster)

#' Project the Adaptive Component Turnover Across the Landscape
#'
#' This function projects the adaptive component turnover across the landscape using Redundancy Analysis (RDA).
#'
#' @param RDA An RDA model object.
#' @param K The number of RDA axes to use for projection.
#' @param env_pres A Raster* or SpatRaster object of environmental variables to project over.
#' @param range An optional RasterLayer or SpatRaster to mask the projections. Default is NULL.
#' @param method The method to use for projection. Either "loadings" or "predict". Default is "loadings".
#' @param scale_env A named vector of scaling factors for the environmental variables.
#' @param center_env A named vector of centering values for the environmental variables.
#'
#' @return A list of RasterLayers representing the projections for the current climates for each RDA axis.
#' @export
#'
adaptive_index <- function(RDA, K, env_pres, range = NULL, method = "loadings", scale_env, center_env){
  
  # Formatting environmental rasters for projection
  var_env_proj_pres <- as.data.frame(rasterToPoints(env_pres[[row.names(RDA$CCA$biplot)]]))
  
  # Standardization of the environmental variables
  var_env_proj_RDA <- as.data.frame(scale(var_env_proj_pres[,-c(1,2)], center_env[row.names(RDA$CCA$biplot)], scale_env[row.names(RDA$CCA$biplot)]))
  
  # Predicting pixels genetic component based on RDA axes
  Proj_pres <- list()
  if(method == "loadings"){
    for(i in 1:K){
      ras_pres <- rasterFromXYZ(data.frame(var_env_proj_pres[,c(1,2)], Z = as.vector(apply(var_env_proj_RDA[,names(RDA$CCA$biplot[,i])], 1, function(x) sum( x * RDA$CCA$biplot[,i])))), crs = crs(env_pres))
      names(ras_pres) <- paste0("RDA_pres_", as.character(i))
      Proj_pres[[i]] <- ras_pres
      names(Proj_pres)[i] <- paste0("RDA", as.character(i))
    }
  }
  
  # Prediction with RDA model and linear combinations
  if(method == "predict"){ 
    pred <- predict(RDA, var_env_proj_RDA[,names(RDA$CCA$biplot[,i])], type = "lc")
    for(i in 1:K){
      ras_pres <- rasterFromXYZ(data.frame(var_env_proj_pres[,c(1,2)], Z = as.vector(pred[,i])), crs = crs(env_pres))
      names(ras_pres) <- paste0("RDA_pres_", as.character(i))
      Proj_pres[[i]] <- ras_pres
      names(Proj_pres)[i] <- paste0("RDA", as.character(i))
    }
  }
  
  # Mask with the range if supplied
  if(!is.null(range)){
    Proj_pres <- lapply(Proj_pres, function(x) mask(x, range))
  }
  
  # Returning projections for current climates for each RDA axis
  return(Proj_pres = Proj_pres)
}

gen <- simple_impute(vcf_to_dosage(liz_vcf))
env <- data.frame(scale(extract(CA_env, liz_coords)))
mod <- rda_run(dos_imputed, env, liz_coords, correctGEO = TRUE)
rda_sig_p <- rda_getoutliers(mod, naxes = "all", outlier_method = "p", p_adj = "fdr", sig = 0.01, plot = FALSE)
rda_snps <- rda_sig_p$rda_snps

## Adaptively enriched RDA (RDA with just outliers)
# ADD GEO CORRECTIPON
RDA_outliers <- vegan::rda(gen[,rda_snps] ~ .,  env)

# The scores of the environmental variables along the RDA axes can be used to calculate a genetic-based index of adaptation for each environmental pixel of the landscape. This index is estimated independently for each RDA axis of interest using the formula:
#$$
#\sum_{i = 1}^{n}{a_ib_i}
#$$
#Where _a_ is the climatic variable score (loading) along the RDA axis, _b_ is the standardized value for this particular variable at the #focal pixel, and _i_ refers to one of the _n_ different variables used in the RDA model.   

## Function to predict the adaptive index across the landscape
## Standardization of the variables
env_scaled <- scale(env, center=TRUE, scale=TRUE) # center=TRUE, scale=TRUE are the defaults for scale()
## Recovering scaling coefficients
scale_env <- attr(env_scaled, 'scaled:scale')
center_env <- attr(env_scaled, 'scaled:center')

range = NULL
res_RDA_proj_current <- adaptive_index(RDA = RDA_outliers, K = 2, env_pres = CA_env, method = "loadings", range = range, scale_env = scale_env, center_env = center_env)

# The adaptive index thus provides an estimate of adaptive genetic similarity or difference of all pixels on the landscape as a function of the values of the environmental predictors at that location. When projected on a map it allows visualizing the different adaptive gradients across a species range.
## Tidying of the climatic rasters for ggplot
res_RDA_proj_current <- rast(stack(res_RDA_proj_current))
RDA_proj <- terra::as.data.frame(res_RDA_proj_current, xy = TRUE)
# For each RDA axes, the adaptive index is normalized between 0 and 1
# Function to normalize a column
normalize <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}
# Apply normalization to all columns except 'x' and 'y'
RDA_proj <- 
  RDA_proj %>%
  dplyr::mutate_at(dplyr::vars(-x, -y), normalize)

## Adaptive genetic turnover projected across lodgepole pine range for RDA1 and RDA2 indexes
TAB_RDA <- as.data.frame(do.call(rbind, RDA_proj[1:2]))
colnames(TAB_RDA)[3] <- "value"
TAB_RDA$variable <- factor(c(rep("RDA1", nrow(RDA_proj[[1]])), rep("RDA2", nrow(RDA_proj[[2]]))), levels = c("RDA1","RDA2"))

TAB_RDA <- RDA_proj %>% pivot_longer(c(-x, -y), names_to = "variable", values_to = "value")

png(here("analysis", "gea", "rda.png"), width = 10, height = 5, units = "in", res = 300)
ggplot(data = TAB_RDA) + 
  geom_raster(aes(x = x, y = y, fill = cut(value, breaks=seq(0, 1, length.out=10), include.lowest = T))) + 
  scale_fill_viridis_d(alpha = 0.8, direction = -1, option = "A", labels = c("Negative scores","","","","Intermediate scores","","","","Positive scores")) +
  xlab("Longitude") + ylab("Latitude") +
  guides(fill=guide_legend(title="Adaptive index")) +
  facet_grid(~ variable) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size=11))
dev.off()

# CLIMATE CHANGE PROJECTION

#### Function to predict genomic offset from a RDA model
genomic_offset <- function(RDA, K, env_pres, env_fut, range = NULL, method = "loadings", scale_env, center_env){
  
  # Mask with the range if supplied
  if(!is.null(range)){
    env_pres <- mask(env_pres, range)
    env_fut <- mask(env_fut, range)
  }
  
  # Formatting and scaling environmental rasters for projection
  var_env_proj_pres <- as.data.frame(scale(rasterToPoints(env_pres[[row.names(RDA$CCA$biplot)]])[,-c(1,2)], center_env[row.names(RDA$CCA$biplot)], scale_env[row.names(RDA$CCA$biplot)]))
  var_env_proj_fut <- as.data.frame(scale(rasterToPoints(env_fut[[row.names(RDA$CCA$biplot)]])[,-c(1,2)], center_env[row.names(RDA$CCA$biplot)], scale_env[row.names(RDA$CCA$biplot)]))

  # Predicting pixels genetic component based on the loadings of the variables
  if(method == "loadings"){
    # Projection for each RDA axis
    Proj_pres <- list()
    Proj_fut <- list()
    Proj_offset <- list()
    for(i in 1:K){
      # Current climates
      ras_pres <- env_pres[[1]]
      ras_pres[!is.na(ras_pres)] <- as.vector(apply(var_env_proj_pres[,names(RDA$CCA$biplot[,i])], 1, function(x) sum( x * RDA$CCA$biplot[,i])))
      names(ras_pres) <- paste0("RDA_pres_", as.character(i))
      Proj_pres[[i]] <- ras_pres
      names(Proj_pres)[i] <- paste0("RDA", as.character(i))
      # Future climates
      ras_fut <- env_fut[[1]]
      ras_fut[!is.na(ras_fut)] <- as.vector(apply(var_env_proj_fut[,names(RDA$CCA$biplot[,i])], 1, function(x) sum( x * RDA$CCA$biplot[,i])))
      Proj_fut[[i]] <- ras_fut
      names(ras_fut) <- paste0("RDA_fut_", as.character(i))
      names(Proj_fut)[i] <- paste0("RDA", as.character(i))
      # Single axis genetic offset 
      Proj_offset[[i]] <- abs(Proj_pres[[i]] - Proj_fut[[i]])
      names(Proj_offset)[i] <- paste0("RDA", as.character(i))
    }
  }
  
  # Predicting pixels genetic component based on predict.RDA
  if(method == "predict"){ 
    # Prediction with the RDA model and both set of envionments 
    pred_pres <- predict(RDA, var_env_proj_pres[,-c(1,2)], type = "lc")
    pred_fut <- predict(RDA, var_env_proj_fut[,-c(1,2)], type = "lc")
    # List format
    Proj_offset <- list()    
    Proj_pres <- list()
    Proj_fut <- list()
    for(i in 1:K){
      # Current climates
      ras_pres <- rasterFromXYZ(data.frame(var_env_proj_pres[,c(1,2)], Z = as.vector(pred_pres[,i])), crs = crs(env_pres))
      names(ras_pres) <- paste0("RDA_pres_", as.character(i))
      Proj_pres[[i]] <- ras_pres
      names(Proj_pres)[i] <- paste0("RDA", as.character(i))
      # Future climates
      ras_fut <- rasterFromXYZ(data.frame(var_env_proj_pres[,c(1,2)], Z = as.vector(pred_fut[,i])), crs = crs(env_pres))
      names(ras_fut) <- paste0("RDA_fut_", as.character(i))
      Proj_fut[[i]] <- ras_fut
      names(Proj_fut)[i] <- paste0("RDA", as.character(i))
      # Single axis genetic offset 
      Proj_offset[[i]] <- abs(Proj_pres[[i]] - Proj_fut[[i]])
      names(Proj_offset)[i] <- paste0("RDA", as.character(i))
    }
  }
  
  # Weights based on axis eigen values
  weights <- RDA$CCA$eig/sum(RDA$CCA$eig)
  
  # Weighing the current and future adaptive indices based on the eigen values of the associated axes
  Proj_offset_pres <- do.call(cbind, lapply(1:K, function(x) rasterToPoints(Proj_pres[[x]])[,-c(1,2)]))
  Proj_offset_pres <- as.data.frame(do.call(cbind, lapply(1:K, function(x) Proj_offset_pres[,x]*weights[x])))
  Proj_offset_fut <- do.call(cbind, lapply(1:K, function(x) rasterToPoints(Proj_fut[[x]])[,-c(1,2)]))
  Proj_offset_fut <- as.data.frame(do.call(cbind, lapply(1:K, function(x) Proj_offset_fut[,x]*weights[x])))
  
  # Predict a global genetic offset, incorporating the K first axes weighted by their eigen values
  ras <- Proj_offset[[1]]
  ras[!is.na(ras)] <- unlist(lapply(1:nrow(Proj_offset_pres), function(x) dist(rbind(Proj_offset_pres[x,], Proj_offset_fut[x,]), method = "euclidean")))
  names(ras) <- "Global_offset"
  Proj_offset_global <- ras
  
  # Return projections for current and future climates for each RDA axis, prediction of genetic offset for each RDA axis and a global genetic offset 
  return(list(Proj_pres = Proj_pres, Proj_fut = Proj_fut, Proj_offset = Proj_offset, Proj_offset_global = Proj_offset_global, weights = weights[1:K]))
}

## Running the function for future climate
env_fut <- CA_env + CA_env^2
res_RDA_future <- genomic_offset(RDA_outliers, K = 2, env_pres = CA_env, env_fut = env_fut, range = range, method = "loadings", scale_env = scale_env, center_env = center_env)

## Tidy data
RDA_proj_offset <- terra::as.data.frame(terra::rast(res_RDA_future$Proj_offset_global), xy = TRUE)

## Projecting genomic offset on a map
colors <- c(
  colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral")[6:5])(2), 
  colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral")[4:3])(2), 
  colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral")[2:1])(3)
)

# Plotting
png(here::here("analysis", "gea", "rda_offset.png"), width = 5, height = 5, units = "in", res = 300)
ggplot(data = RDA_proj_offset) + 
  geom_raster(aes(x = x, y = y, fill = cut(Global_offset, breaks = seq(1, 8, by = 1), include.lowest = TRUE)), alpha = 1) + 
  scale_fill_manual(
    values = colors, 
    labels = c("1-2", "2-3", "3-4", "4-5", "5-6", "6-7", "7-8"), 
    guide = guide_legend(
      title = "Genomic offset", 
      title.position = "top", 
      title.hjust = 0.5, 
      ncol = 1, 
      label.position = "right"
    ), 
    na.translate = FALSE
  ) +
  xlab("Longitude") + ylab("Latitude") +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(
    panel.grid = element_blank(), 
    plot.background = element_blank(), 
    panel.background = element_blank(), 
    strip.text = element_text(size = 11)
  )
dev.off()