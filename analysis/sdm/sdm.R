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
sdm <- function(coordinates, envlayers, biasdat = NULL, models = "MAXENT.Phillips", nbkg = 10000, output = FALSE, file.name = "", tune = TRUE){
  
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
                                       PA.nb.rep = 0) #number of required Pseudo Absences selection (if needed). 0 by Default.
  
  
  if (tune) {
    # Tune parameters
    Biomod.tuning <- BIOMOD_Tuning(myBiomodData,
                                   models = models)
    
    myBiomodOption <- Biomod.tuning$models.options
  } else {
    #default parameters
    myBiomodOption <- BIOMOD_ModelingOptions()
  }
  
  # Setting up for cross evaluation
  DataSplitTable <- BIOMOD_cv(myBiomodData, k = 10, repetition = 1, do.full.models = T)
  
  # SDM model
  myBiomodModelOut <- BIOMOD_Modeling( myBiomodData,
                                       models = models, #choosing model
                                       models.options = myBiomodOption,
                                       NbRunEval = 1, #CHANGE THIS LATER
                                       DataSplitTable = DataSplitTable,
                                       Prevalence = NULL, #used to create weighted responses (NULL=no weighting)
                                       VarImport = 0, #number of evals to run to determine variable importance
                                       models.eval.meth = c('TSS', 'ROC'), #evaluation methods
                                       SaveObj = TRUE,
                                       #if true, all model prediction will be scaled with a binomial GLM
                                       rescal.all.models = FALSE, #model prediction will be scaled with a binomial GLM
                                       do.full.models = TRUE, # (ADVISED AGAINST IN DOCS) create model calibrated and evaluated with the whole dataset
                                       modeling.id = "test")
  
  # Get evaluation statistics for each model
  myBiomodModelEval <- get_evaluations(myBiomodModelOut)
  
  # Display the model evaluation statistics
  myBiomodModelEval["TSS", "Testing.data" ,,,] #True Skill Statistic
  myBiomodModelEval["ROC", "Testing.data",,,] #Same as AUC (Area under Receiver Operating Curve(ROC))
  
  # Mean and sd of results
  eval.results <- data.frame(TSS_mean=mean(myBiomodModelEval["TSS", "Testing.data" ,,,]),
                           TSS_min=min(myBiomodModelEval["TSS", "Testing.data" ,,,]),
                           TSS_sd=sd(myBiomodModelEval["TSS", "Testing.data" ,,,]),
                           AUC_mean=mean(myBiomodModelEval["ROC", "Testing.data",,,]),
                           AUC_min=min(myBiomodModelEval["ROC", "Testing.data" ,,,]),
                           AUC_sd=sd(myBiomodModelEval["ROC", "Testing.data",,,]))
  eval.results <- round(eval.results, 2)
  eval.results
  
  
  # Project model using final model (Run 11)
  myBiomodProj <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
                                    new.env = envlayers, # project
                                    #selected.models = paste0(BIOMOD_name, "_AllData_RUN11_", models), #which models to use
                                    proj.name = 'BIOMOD', # name of projection
                                    binary.meth = NULL, # evaluation method by which to choose threshold (If NULL then no binary transformation computed, else the given binary techniques will be used to transform the projection into 0/1 data.)
                                    compress = 'xz', # compression format for object storage
                                    build.clamping.mask = FALSE, #if TRUE, a clamping mask will be saved on hard drive
                                    output.format = '.grd')

  lh_hindcast <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
                                    new.env = lh, # project
                                    #selected.models = paste0(BIOMOD_name, "_AllData_RUN11_", models), #which models to use
                                    proj.name = 'BIOMOD', # name of projection
                                    binary.meth = NULL, # evaluation method by which to choose threshold (If NULL then no binary transformation computed, else the given binary techniques will be used to transform the projection into 0/1 data.)
                                    compress = 'xz', # compression format for object storage
                                    build.clamping.mask = FALSE, #if TRUE, a clamping mask will be saved on hard drive
                                    output.format = '.grd')  
  lgm_hindcast <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
                                    new.env = lgm, # project
                                    #selected.models = paste0(BIOMOD_name, "_AllData_RUN11_", models), #which models to use
                                    proj.name = 'BIOMOD', # name of projection
                                    binary.meth = NULL, # evaluation method by which to choose threshold (If NULL then no binary transformation computed, else the given binary techniques will be used to transform the projection into 0/1 data.)
                                    compress = 'xz', # compression format for object storage
                                    build.clamping.mask = FALSE, #if TRUE, a clamping mask will be saved on hard drive
                                    output.format = '.grd')                                     
  # Plot SDM
  # COME BACK TO THIS MEAN THING: 
  sdm_raster_BIOMOD <- mean(get_predictions(myBiomodProj))
  lh <- mean(get_predictions(lh_hindcast))/1000
  lgm <- mean(get_predictions(lgm_hindcast))/1000
  
  # Convert raster from 0 to 1000 to 0 to 1 
  # (biomod automatically multiplies everything by 1000 to make integer values)
  sdm_raster <- sdm_raster_BIOMOD/1000 
  par(mfrow = c(1,3))
  plot(sdm_raster, col = viridis::inferno(100), axes = FALSE, box = FALSE)
  points(coordinates, col = "red")
  plot(lh, col = viridis::inferno(100), axes = FALSE, box = FALSE)
  points(coordinates, col = "red")
  plot(lgm, col = viridis::inferno(100), axes = FALSE, box = FALSE)
  points(coordinates, col = "red")

  # Make outputs into list
  spp.eval = list(raster = sdm_raster, eval = eval.results)
  
  # Write out evaluation results
  if (output) write.csv(spp.eval[["eval"]], here("outputs", paste0(file.name, "_eval_results.csv")))
  
  # Delete automatically generated biomod directory
  unlink(BIOMOD_name, recursive = TRUE) 
  
  return(spp.eval)
}

