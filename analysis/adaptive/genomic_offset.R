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
  
  # Deal with future layer naming
  env_fut_1 <- raster::subset(env_fut, c(1,3))
  names(env_fut_1) <- names(env_pres)
  env_fut_2 <- raster::subset(env_fut, 2:3)
  names(env_fut_2) <- names(env_pres)
  
  # Formatting and scaling environmental rasters for projection
  var_env_proj_pres <- as.data.frame(scale(raster::rasterToPoints(env_pres[[row.names(biplot)]])[,-c(1,2)], center_env[row.names(biplot)], scale_env[row.names(biplot)]))
  var_env_proj_fut_1 <- as.data.frame(scale(raster::rasterToPoints(env_fut_1[[row.names(biplot)]])[,-c(1,2)], center_env[row.names(biplot)], scale_env[row.names(biplot)]))
  var_env_proj_fut_2 <- as.data.frame(scale(raster::rasterToPoints(env_fut_2[[row.names(biplot)]])[,-c(1,2)], center_env[row.names(biplot)], scale_env[row.names(biplot)]))
  
  # Predicting pixels genetic component based on the loadings of the variables
  if(method == "loadings"){
    # Projection for each RDA axis
    Proj_pres <- offset_proj_helper(env = env_pres, biplot = biplot, var_env_proj = var_env_proj_pres, K = K, type = "present")
    Proj_fut_1 <- offset_proj_helper(env = env_fut_1, biplot = biplot, var_env_proj = var_env_proj_fut_1, K = K, type = "future")
    Proj_fut_2 <- offset_proj_helper(env = env_fut_2, biplot = biplot, var_env_proj = var_env_proj_fut_2, K = K, type = "future")

    Proj_offset_1 <- list()
    Proj_offset_2 <- list()
    
    # Single axis genetic offset 
    for(i in 1:K) {
      Proj_offset_1[[i]] <- abs(Proj_pres[[i]] - Proj_fut_1[[i]])
      names(Proj_offset_1)[i] <- paste0("RDA", as.character(i))
      Proj_offset_2[[i]] <- abs(Proj_pres[[i]] - Proj_fut_2[[i]])
      names(Proj_offset_2)[i] <- paste0("RDA", as.character(i))
    }
  }

  # # Predicting pixels genetic component based on predict.RDA
  # if (method == "predict") { 
  #   # Prediction with the RDA model and both set of environments 
  #   pred_pres <- predict(mod, var_env_proj_pres[,-c(1,2)], type = "lc")
  #   pred_fut <- predict(mod, var_env_proj_fut[,-c(1,2)], type = "lc")
  #   # List format
  #   Proj_offset <- list()    
  #   Proj_pres <- list()
  #   Proj_fut <- list()
  #   for(i in 1:K) {
  #     # Current climates
  #     ras_pres <- rasterFromXYZ(data.frame(var_env_proj_pres[,c(1,2)], Z = as.vector(pred_pres[,i])), crs = crs(env_pres))
  #     names(ras_pres) <- paste0("RDA_pres_", as.character(i))
  #     Proj_pres[[i]] <- ras_pres
  #     names(Proj_pres)[i] <- paste0("RDA", as.character(i))
  #     # Future climates
  #     ras_fut <- rasterFromXYZ(data.frame(var_env_proj_pres[,c(1,2)], Z = as.vector(pred_fut[,i])), crs = crs(env_pres))
  #     names(ras_fut) <- paste0("RDA_fut_", as.character(i))
  #     Proj_fut[[i]] <- ras_fut
  #     names(Proj_fut)[i] <- paste0("RDA", as.character(i))
  #     # Single axis genetic offset 
  #     Proj_offset[[i]] <- abs(Proj_pres[[i]] - Proj_fut[[i]])
  #     names(Proj_offset)[i] <- paste0("RDA", as.character(i))
  #   }
  # }

  # Weights based on axis eigenvalues
  weights <- eig %>% dplyr::mutate(weights = mod.CCA.eig / sum(eig$mod.CCA.eig)) %>% pull(weights)

  # Weighing the current and future adaptive indices based on the eigenvalues of the associated axes
  Proj_offset_pres <- do.call(cbind, lapply(1:K, function(x) rasterToPoints(Proj_pres[[x]])[,-c(1,2)]))
  Proj_offset_pres <- as.data.frame(do.call(cbind, lapply(1:K, function(x) Proj_offset_pres[,x] * weights[x])))

  Proj_offset_fut_1 <- do.call(cbind, lapply(1:K, function(x) rasterToPoints(Proj_fut_1[[x]])[,-c(1,2)]))
  Proj_offset_fut_1 <- as.data.frame(do.call(cbind, lapply(1:K, function(x) Proj_offset_fut_1[,x] * weights[x])))

  Proj_offset_fut_2 <- do.call(cbind, lapply(1:K, function(x) rasterToPoints(Proj_fut_2[[x]])[,-c(1,2)]))
  Proj_offset_fut_2 <- as.data.frame(do.call(cbind, lapply(1:K, function(x) Proj_offset_fut_2[,x] * weights[x])))
  
  # Predict a global genetic offset, incorporating the K first axes weighted by their eigenvalues
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
