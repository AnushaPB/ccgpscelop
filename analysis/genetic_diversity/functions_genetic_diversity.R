

get_het <- function(file = "58-Sceloporus_chr.het"){
  # callable sites for denominator
  callable <- read.csv(here("data_processing", "callable_counts.csv"))
  callable_sites <- callable %>% pull(callable_sites)

  het <- format_het(here("analysis", "genetic_diversity", "outputs", file), callable_sites = callable_sites)
  
  if ("IID" %in% names(het)) het$SampleID <- het$IID
  if ("INDV" %in% names(het)) het$SampleID <- het$INDV

  return(het)
}

format_het <- function(path, callable_sites){
  # Load the data in R
  het_data <- read_table(path)

  # Calculate the average heterozygosity per individual
  if ("N(NM)" %in% names(het_data)) het_data$Ho <- (het_data$`N(NM)` - het_data$`O(HOM)`)/callable_sites
  if ("N_SITES" %in% names(het_data)) het_data$Ho <- (het_data$N_SITES - het_data$`O(HOM)`)/callable_sites
  return(het_data) 
}

get_roh <- function(){
  coords <- get_coords(sf = TRUE)
  pops <- get_pops()
  coords <- left_join(coords, pops, by = "SampleID")

  roh <-
    map(paste0("pop", 1:5), ~{
      poproh <- 
        read_table(here("analysis", "genetic_diversity", "outputs", paste0("K5_", .x, ".froh")), col_names = FALSE) %>%
        rename(SampleID = X1, froh = X2) %>%
        right_join(filter(coords, pop == .x)) %>%
        # missing froh values indicates no roh greater than the minimum size were found
        # remains as NA if pop was not found (i.e., individual was not included in analysis)
        mutate(
          froh0 = case_when(is.na(froh) & !is.na(pop) ~ 0, .default = froh),
          pop = .x
        )
    }) %>% bind_rows()

  return(roh)
}

# Impute NA values for extracted values using bilinear interpolation
# x - imput vector of extracted values
# r - raster from which values were extracted
# coords - coordinates of the points where values were extracted
# Checks which values are NA and replaces them with bilinear interpolation from the raster
bilinear_impute <- function(x, r, coords) {
  na_vals <- which(is.na(x))
  if (length(na_vals) > 0) {
    message("Imputing ", length(na_vals), " missing values")
    x[na_vals] <- terra::extract(r, coords[na_vals,], ID = FALSE, method = "bilinear")[,1]
  }
  return(x)
}

geosummarize <- function(coords, stat = "Ho", res = 50000){
  # transform coords
  coords <- coords %>% st_transform(3310)

  # Load or create a raster
  lyr <- wingen::coords_to_raster(coords, res = res)

  # Get raster coords
  raster_coords <- terra::extract(lyr, coords, xy = TRUE)

  # combine with coords
  final_coords <- 
    bind_cols(raster_coords, coords) %>%
    group_by(x, y) %>%
    summarize_at({{stat}}, mean, na.rm = TRUE) %>%
    drop_na(x, y) %>%
    st_as_sf(coords = c("x", "y"), crs = st_crs(coords)) %>%
    st_transform(4326)

}

spatial_dredge <- function(full_formula, data, coords = c("x", "y"), 
                          random = "~ 1 | dummy", method = "ML", corFunction = corExp) {
  
  # Get all predictor variables from formula
  predictors <- all.vars(full_formula[[3]])
  
  # Create all possible combinations of predictors
  n_preds <- length(predictors)
  combinations <- lapply(1:n_preds, 
                        function(k) combn(predictors, k, simplify = FALSE))
  combinations <- unlist(combinations, recursive = FALSE)

  # Add null model (intercept only)
  combinations <- c(list(character(0)), combinations)
  
  # Function to fit a model and return its metrics
  fit_model <- function(predictors) {
    tryCatch({
      # Construct formula
      if(length(predictors) == 0) {
        # Null model (intercept only)
        current_formula <- as.formula(paste("Ho ~ 1"))
      } else {
        # Model with predictors
        current_formula <- as.formula(
          paste("Ho ~", paste(predictors, collapse = " + "))
        )
      }
      
      # Fit spatial model
      current_model <- lme(fixed = current_formula,
                          data = data,
                          random = as.formula(random),
                          correlation = corFunction(1, form = as.formula(
                            paste("~", paste(coords, collapse = " + "))
                          ), nugget = TRUE),
                          method = method)
      
      # Return results
      return(data.frame(
        model_formula = paste(deparse(current_formula), collapse = ""),
        AIC = AIC(current_model),
        BIC = BIC(current_model),
        logLik = logLik(current_model),
        stringsAsFactors = FALSE
      ))
      
    }, error = function(e) {
      cat("Error fitting model:", deparse(current_formula), "\n")
      cat("Error message:", conditionMessage(e), "\n")
      return(NULL)
    })
  }
  
  # Run models for each combination using purrr::map
  model_results <- future_map_dfr(combinations, fit_model, .progress = TRUE)

  # Sort results by AIC
  model_results <- 
    model_results %>%
    arrange(AIC)
  
  return(model_results)
}

make_listw <- function(coords, nbdist = 10000){
  # Create a spatial weights matrix
  listw <- 
    spdep::nb2listw(
      spdep::dnearneigh(
        coords, 
        d1 = 0, 
        d2 = nbdist, 
        longlat = TRUE
      ),
      zero.policy = TRUE
    )

  return(listw)
  
}

sem_dredge <- function(full_formula, data, listw) {
  
  # Get all predictor variables from formula
  predictors <- all.vars(full_formula[[3]])
  
  # Create all possible combinations of predictors
  n_preds <- length(predictors)
  combinations <- lapply(1:n_preds, 
                        function(k) combn(predictors, k, simplify = FALSE))
  combinations <- unlist(combinations, recursive = FALSE)

  # Add null model (intercept only)
  combinations <- c(list(character(0)), combinations)
  
  # Function to fit a model and return its metrics
  fit_model <- function(predictors) {
    tryCatch({
      # Construct formula
      if(length(predictors) == 0) {
        # Null model (intercept only)
        current_formula <- as.formula(paste("Ho ~ 1"))
      } else {
        # Model with predictors
        current_formula <- as.formula(
          paste("Ho ~", paste(predictors, collapse = " + "))
        )
      }
      
      # Fit spatial model
      current_model <- errorsarlm(current_formula, data = data, listw = listw, na.action = "na.omit", zero.policy = TRUE)
      
      # Return results
      return(data.frame(
        model_formula = paste(deparse(current_formula), collapse = ""),
        AIC = AIC(current_model),
        BIC = BIC(current_model),
        logLik = logLik(current_model),
        stringsAsFactors = FALSE
      ))
      
    }, error = function(e) {
      cat("Error fitting model:", deparse(current_formula), "\n")
      cat("Error message:", conditionMessage(e), "\n")
      return(NULL)
    })
  }
  
  # Run models for each combination using purrr::map
  model_results <- future_map_dfr(combinations, fit_model, .progress = TRUE)

  # Sort results by AIC
  model_results <- 
    model_results %>%
    arrange(AIC)
  
  return(model_results)
}

lm_variance_partition <- function(predictor_sets, response, scale_vars = TRUE) {
  # Validate inputs
  if (length(predictor_sets) < 2) {
    stop("At least two predictor sets are required")
  }
  
  # Convert everything to data frames
  predictor_sets <- map(predictor_sets, as.data.frame)
  
  # Validate dimensions
  n_obs <- length(response)
  walk2(
    predictor_sets,
    names(predictor_sets),
    ~if (nrow(.x) != n_obs) {
      stop(sprintf("Predictor set %s has different number of observations than response", .y))
    }
  )
  
  # Scale variables if requested
  if (scale_vars) {
    response <- scale(response) %>% as.vector()
    predictor_sets <- map(predictor_sets, ~scale(.x) %>% as.data.frame())
  }
  
  # Function to calculate R2 and adjusted R2
  get_r2 <- function(model) {
    sum <- summary(model)
    c(
      r2 = sum$r.squared,
      adj_r2 = sum$adj.r.squared
    )
  }
  
  # Individual models for each predictor set
  individual_models <- predictor_sets %>%
    map(~lm(response ~ ., data = .x))
  
  # Full model with all predictors
  all_predictors <- as.data.frame(predictor_sets)
  names(all_predictors) <- names(predictor_sets)
  full_model <- lm(response ~ ., data = all_predictors)
  
  # Get R2 values for individual models
  individual_r2 <- map(individual_models, get_r2) %>%
    bind_rows(.id = "predictor_set")
  
  # Function to get all possible combinations
  get_combinations <- function(n) {
    combn(names(predictor_sets), n, simplify = FALSE)
  }
  
  # Get all combinations of predictor sets
  all_combinations <- map(
    2:length(predictor_sets),
    get_combinations
  ) %>%
    flatten()
  
  # Calculate R2 for all combinations
  combination_results <- all_combinations %>%
    set_names(map_chr(., ~paste(.x, collapse = "+"))) %>%
    map(~{
      combined_predictors <- predictor_sets[.x] %>% as.data.frame()
      names(combined_predictors) <- .x
      model <- lm(response ~ ., data = combined_predictors)
      get_r2(model)
    }) %>%
    bind_rows(.id = "combination")
  
  # Calculate unique contributions (Type III SS)
  unique_contributions <- map2_dfr(
    names(predictor_sets),
    predictor_sets,
    ~{
      # Model without this predictor set
      other_predictors <- predictor_sets[names(predictor_sets) != .x] %>% as.data.frame()
      names(other_predictors) <- names(predictor_sets)[names(predictor_sets) != .x]
      reduced_model <- lm(response ~ ., data = other_predictors)
      
      # Calculate unique contribution
      unique_r2 <- get_r2(full_model)[2] - get_r2(reduced_model)[2]
      
      tibble(
        predictor_set = .x,
        unique_contribution = unique_r2
      )
    }
  )
  
  # Calculate shared variance
  total_r2 <- get_r2(full_model)[2]
  shared_variance <- total_r2 - sum(unique_contributions$unique_contribution)
  
  # Prepare summary statistics
  summary_stats <- tibble(
    total_explained = total_r2,
    shared = shared_variance,
    unexplained = 1 - total_r2
  )
  
  
  # Print summary
  cat("\nVariance Partitioning Summary:\n")
  cat("----------------------------\n")
  ind <- 
    individual_r2 %>% 
    mutate(across(where(is.numeric), ~scales::percent(.x, accuracy = 1)))
  uniq <- 
    unique_contributions %>% 
    mutate(unique_contribution = scales::percent(unique_contribution, accuracy = 1))
  summary_table <- left_join(ind, uniq, by = "predictor_set")
  print(summary_table)

  cat("\nSummary statistics:\n")
  print(summary_stats %>% 
        mutate(across(everything(), ~scales::percent(.x, accuracy = 1))))
  
  # Create results list
  results <- list(
    summary_table = summary_table,
    individual_models = individual_models,
    full_model = full_model,
    individual_r2 = individual_r2,
    unique_contributions = unique_contributions,
    combination_results = combination_results,
    summary_stats = summary_stats
  )
  return(results)
}

