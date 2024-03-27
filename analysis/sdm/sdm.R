
#' Create an SDM using BIOMOD
#'
#' @param coordinates coordinates of species occurrence
#' @param envlayers environmental layers to use for modelling
#' @param biasdat optional bias data (coordinates)
#' @param models SDM models to run 
#' @param nbkg number of background points to use if `biasdat` is not provided
#' @param output whether to save evaluation output as a csv
#' @param file.name file name for evaluation output
#' @param tune whether to tune the model parameters (otherwise uses the defaults)
#'
#' @return
#' @export
#'
#' @examples
sdm <- function(coordinates, envlayers, biasdat = NULL, models = "MAXENT", nbkg = 10000, output = FALSE, file.name = "", tune = TRUE){
  
  # Create a name for BIOMOD files (using random string of numbers so that each is unique if running in parallel)
  BIOMOD_name <- paste0("BIOMOD", paste(sample(0:100, 5, replace = TRUE), collapse = ""))

  # Create random background points or use bias data if provided
  if(is.null(biasdat)){
    # Generate random background points
    abs <- sampleRandom(envlayers[[1]], nbkg, xy = TRUE)
    bk.df <- data.frame(pa = 0, abs[, c("x","y")])
  } else {
    # Create presence/absence data frame using bias data
    bk.df <- data.frame(pa = 0, biasdat)
    colnames(bk.df) <- c("pa", "x", "y")
  }
  
  # Create presence/absence dataframe
  pa.df <- data.frame(pa = 1, coordinates)
  colnames(pa.df) <- c("pa", "x", "y")
  pa.df <- rbind(pa.df, bk.df)
  
  # Format Resp
  # Response variable (presence/absence vector)
  myResp <- as.numeric(pa.df[,"pa"]) 
  # Coordinates for response variable
  myRespXY <- pa.df[,c("x","y")] 
  
  # Create object (myBiomodData) to contain all the previous objects within it,formatted correctly
  myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,#presence absence vector
                                       expl.var = envlayers, #environmental layers
                                       resp.xy = myRespXY, #coordinates                                  
                                       resp.name = BIOMOD_name, # name for files - doesn't really matter 
                                       PA.nb.rep = 0, #number of required Pseudo Absences selection (if needed). 0 by Default.
                                       ) 
  
  if (tune) {
    # Tune parameters
    Biomod.tuning <- BIOMOD_Tuning(myBiomodData,
                                   models = models)
    
    myBiomodOptions <- Biomod.tuning$models.options
  } else {
    #default parameters
    myBiomodOptions <- BIOMOD_ModelingOptions()
  }
  
  # SDM model
  myBiomodModelOut <- BIOMOD_Modeling(bm.format = myBiomodData,
                                    modeling.id = BIOMOD_name,
                                    models = models,
                                    bm.options = myBiomodOptions,
                                    CV.strategy = 'kfold',
                                    CV.nb.rep = 1,
                                    CV.k = 10,
                                    CV.do.full.models = TRUE,
                                    metric.eval = c('TSS','ROC'),
                                    #var.import = 2, # number of permutations for variable importance
                                    seed.val = 42)
  
  # Get evaluation statistics for each model
  # Contains evaluation metric for different models and dataset. 
  # Evaluation metric are calculated on the calibrating data (column calibration), 
  # on the cross-validation data (column validation) or on the 
  # evaluation data (column evaluation). We care about the validation data
  myBiomodModelEval <- get_evaluations(myBiomodModelOut)
  
  # Mean
  eval.results <- 
    myBiomodModelEval %>% 
    dplyr::group_by(metric.eval) %>% 
    dplyr::summarize_at(c("sensitivity", "specificity", "validation"), mean, na.rm = TRUE)
  print(eval.results)
  
  # Project model using final model (Run 11)
  myBiomodProj <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
                                    new.env = envlayers, # project
                                    #selected.models = paste0(BIOMOD_name, "_AllData_RUN11_", models), #which models to use
                                    proj.name = 'BIOMOD', # name of projection
                                    binary.meth = NULL, # evaluation method by which to choose threshold (If NULL then no binary transformation computed, else the given binary techniques will be used to transform the projection into 0/1 data.)
                                    compress = 'xz', # compression format for object storage
                                    build.clamping.mask = FALSE, #if TRUE, a clamping mask will be saved on hard drive
                                    output.format = '.grd')                                
  # Get predictions from full model 
  predictions <- get_predictions(myBiomodProj)
  allRun <- grepl("allData_allRun", names(predictions))
  sdm_raster_BIOMOD <- predictions[[allRun]]
  
  # Convert raster from 0 to 1000 to 0 to 1 
  # (biomod automatically multiplies everything by 1000 to make integer values)
  sdm_raster <- sdm_raster_BIOMOD/1000 

  # Make outputs into list
  spp.eval = list(raster = sdm_raster, eval = eval.results)
  
  # Write out evaluation results
  if (output) write.csv(spp.eval[["eval"]], here("outputs", paste0(file.name, "_eval_results.csv")))
  
  # Delete automatically generated biomod directory
  # DON"T DO THIS IF YOU WANT TO HINDCAST
  #unlink(BIOMOD_name, recursive = TRUE) 
  
  return(list(eval = spp.eval, myBiomodModelOut = myBiomodModelOut, raster = sdm_raster))
}

hindcast <- function(myBiomodModelOut, new_env){
  hindcast_proj <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
                                      new.env = new_env, # project
                                      proj.name = 'BIOMOD', # name of projection
                                      binary.meth = NULL, # evaluation method by which to choose threshold (If NULL then no binary transformation computed, else the given binary techniques will be used to transform the projection into 0/1 data.)
                                      compress = 'xz', # compression format for object storage
                                      build.clamping.mask = FALSE, #if TRUE, a clamping mask will be saved on hard drive
                                      output.format = '.grd')                                    
  # Plot SDM
  hindcast_raster <- get_predictions(hindcast_proj)
  allRun <- grepl("allData_allRun", names(hindcast_raster))
  # Convert raster from 0 to 1000 to 0 to 1 
  # (biomod automatically multiplies everything by 1000 to make integer values)
  hindcast_raster <- hindcast_raster[[allRun]]/1000

  return(hindcast_raster)
}

pc_var_selection <- function(envlayers, coords = NULL){
  # If coordinates are not provided, get them
  if (is.null(coords)) coords <- get_coords(sf = TRUE)

  # Perform Principal Component Analysis (PCA) on environmental layers
  pc <- rasterPCA(envlayers)

  # Get loadings from PCA model
  loadings <- pc$model$loadings

  # Get indices of top variables based on maximum absolute loading values
  top_vars_indices <- apply(loadings, 2, function(x) which.max(abs(x)))

  # Get names of top variables
  top_vars_names <- rownames(loadings)[top_vars_indices]

  # Extract values of top variables at given coordinates
  cur_vals <- extract(envlayers, coords)

  # Remove duplicate names if present, assuming order matters and first occurrence is priority
  top_vars_names <- unique(top_vars_names)

  # Initialize list of selected variables
  selected_vars <- c()

  for (var in top_vars_names) {
    # Skip if the variable isn't in the dataframe
    if (!var %in% names(data)) next
    
    # If it's the first variable, just add it
    if (length(selected_vars) == 0) {
      selected_vars <- c(selected_vars, var)
    } else {
      # Check correlation with already selected variables
      cor_values <- sapply(selected_vars, function(x) abs(cor(cur_vals[,x], cur_vals[,var], use = "complete.obs")))
      
      # Add the variable if it's not highly correlated with any of the selected variables
      if (all(cor_values <= 0.6)) {
        selected_vars <- c(selected_vars, var)
      }
    }
  }

  return(selected_vars)
}
