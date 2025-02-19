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
  
  # Formatting and scaling environmental rasters for projection
  var_env_proj_pres <- as.data.frame(scale(raster::rasterToPoints(env_pres[[row.names(biplot)]])[,-c(1,2)], center_env[row.names(biplot)], scale_env[row.names(biplot)]))
  var_env_proj_fut <- as.data.frame(scale(raster::rasterToPoints(env_fut[[row.names(biplot)]])[,-c(1,2)], center_env[row.names(biplot)], scale_env[row.names(biplot)]))
  
  # Predicting pixels genetic component based on the loadings of the variables
  if(method == "loadings"){
    # Projection for each RDA axis
    Proj_pres <- list()
    Proj_fut <- list()
    Proj_offset <- list()
    
    for(i in 1:K) {
      # Current climates
      ras_pres <- env_pres[[1]]
      ras_pres[!is.na(ras_pres)] <- as.vector(apply(var_env_proj_pres[,rownames(biplot[i])], 1, function(x) sum(x * biplot[,i])))
      names(ras_pres) <- paste0("RDA_pres_", as.character(i))
      Proj_pres[[i]] <- ras_pres
      names(Proj_pres)[i] <- paste0("RDA", as.character(i))
      
      # Future climates
      ras_fut <- env_fut[[1]]
      ras_fut[!is.na(ras_fut)] <- as.vector(apply(var_env_proj_fut[,rownames(biplot[i])], 1, function(x) sum(x * biplot[,i])))
      Proj_fut[[i]] <- ras_fut
      names(ras_fut) <- paste0("RDA_fut_", as.character(i))
      names(Proj_fut)[i] <- paste0("RDA", as.character(i))
      
      # Single axis genetic offset 
      Proj_offset[[i]] <- abs(Proj_pres[[i]] - Proj_fut[[i]])
      names(Proj_offset)[i] <- paste0("RDA", as.character(i))
    }
  }
  
  # Predicting pixels genetic component based on predict.RDA
  if (method == "predict") { 
    # Prediction with the RDA model and both set of environments 
    pred_pres <- predict(mod, var_env_proj_pres[,-c(1,2)], type = "lc")
    pred_fut <- predict(mod, var_env_proj_fut[,-c(1,2)], type = "lc")
    # List format
    Proj_offset <- list()    
    Proj_pres <- list()
    Proj_fut <- list()
    for(i in 1:K) {
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

  # Weights based on axis eigenvalues
  weights <- eig %>% dplyr::mutate(weights = mod.CCA.eig / sum(eig$mod.CCA.eig)) %>% pull(weights)

  # Weighing the current and future adaptive indices based on the eigenvalues of the associated axes
  Proj_offset_pres <- do.call(cbind, lapply(1:K, function(x) rasterToPoints(Proj_pres[[x]])[,-c(1,2)]))
  Proj_offset_pres <- as.data.frame(do.call(cbind, lapply(1:K, function(x) Proj_offset_pres[,x] * weights[x])))
  Proj_offset_fut <- do.call(cbind, lapply(1:K, function(x) rasterToPoints(Proj_fut[[x]])[,-c(1,2)]))
  Proj_offset_fut <- as.data.frame(do.call(cbind, lapply(1:K, function(x) Proj_offset_fut[,x] * weights[x])))
  
  # Predict a global genetic offset, incorporating the K first axes weighted by their eigenvalues
  ras <- Proj_offset[[1]]
  ras[!is.na(ras)] <- unlist(lapply(1:nrow(Proj_offset_pres), function(x) dist(rbind(Proj_offset_pres[x,], Proj_offset_fut[x,]), method = "euclidean")))
  names(ras) <- "Global_offset"
  Proj_offset_global <- ras
  
  # Return projections for current and future climates for each RDA axis, prediction of genetic offset for each RDA axis and a global genetic offset 
  return(list(Proj_pres = Proj_pres, Proj_fut = Proj_fut, Proj_offset = Proj_offset, Proj_offset_global = Proj_offset_global, weights = weights[1:K]))
}