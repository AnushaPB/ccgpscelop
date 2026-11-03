suppressMessages({
  library(tidyverse)
  library(raster)
  # library(terra)
  library(vegan)
  library(qvalue)
  library(peakRAM)
  library(here)
  library(devtools)
  library(algatr)
})

# if (!require("algatr", character.only = TRUE)) {
#   # Install the package if not installed
#   devtools::install_github("TheWangLab/algatr", quiet = T)
# }

#set up log file writing
log_smk <- function() {
  if (exists("snakemake") & length(snakemake@log) != 0) {
    log <- file(snakemake@log[1][[1]], open = "wt")
    sink(log, append = TRUE)
    sink(log, append = TRUE, type = "message")
  }
}

log_smk()

#!/usr/bin/env Rscript
species = snakemake@params[[1]]
output_path = snakemake@params[[2]]
rmislands = as.logical(snakemake@params[[3]])
impute = snakemake@params[[4]]
correctGEO = as.logical(snakemake@params[[5]])
correctPC = as.logical(snakemake@params[[6]]) # true or false
pc_path = snakemake@params[[7]]
nPC = snakemake@params[[8]]
pc_data_type = snakemake@params[[9]]
species_ks = snakemake@params[[10]]
model = snakemake@params[[11]]
sig = snakemake@params[[12]] # only for outlier_method = "p"
z = snakemake@params[[13]]# only for outlier_method = "z"
Pin = snakemake@params[[14]]
R2permutations = snakemake@params[[15]]
R2scope = as.logical(snakemake@params[[16]])
p_adj = snakemake@params[[17]] # only if outlier_method = "p"
save_impute = as.logical(snakemake@params[[18]])
intervals = as.logical(snakemake@params[[19]])
scaff = as.character(snakemake@params[[20]])
env_var_type = snakemake@params[[21]]
shape_path = as.character(snakemake@params[[22]])

vcf = snakemake@input[["vcf"]]

source(paste0(snakemake@scriptdir, "/general_functions.R"))


# Import and process data -------------------------------------------------
peakRAM_imp <-
  peakRAM::peakRAM(
    dat <- get_input_objects(species = species, 
                             analysis = "vcf",
                             impute = impute,
                             rmislands = rmislands,
                             save_impute = FALSE, # save this for later
                             intervals = intervals,
                             scaff = scaff,
                             incl_env = TRUE,
                             env_var_type = env_var_type,
                             vcf_path = vcf,
                             shape_path = shape_path)
  )

# Explicit outputs --------------------------------------------------------

zscore_output = snakemake@output[["zscore_output"]]
rdadapt_output = snakemake@output[["rdadapt_output"]]
cortest_output = snakemake@output[["cortest_output"]]
imputed_output = snakemake@output[["imputed_output"]]
peakram_output = snakemake@output[["peakram_output"]]
colsum_output = snakemake@output[["colsum_output"]]
ybar_output = snakemake@output[["ybar_output"]]
v_output = snakemake@output[["v_output"]]
u_output = snakemake@output[["u_output"]]
wa_output = snakemake@output[["wa_output"]]
qr_output = snakemake@output[["qr_output"]]
eig_output = snakemake@output[["eig_output"]]
biplot_output = snakemake@output[["biplot_output"]]
qraux_output = snakemake@output[["qraux_output"]]
envcentre_output = snakemake@output[["envcentre_output"]]
chi_output = snakemake@output[["chi_output"]]
scalload_output = snakemake@output[["scalload_output"]]
unscalload_output = snakemake@output[["unscalload_output"]]


# If no SNPs in scaffold, touch output files and stop ---------------------

if (is.null(ncol(dat$gen))) {
  # Create all output files
  file.create(# Model-related output files
    colsum_output, ybar_output, v_output, u_output, wa_output, qr_output,
    eig_output, biplot_output, qraux_output, envcentre_output, chi_output,
    scalload_output, unscalload_output,
    # Outlier-related output files
    zscore_output, rdadapt_output,
    # Correlation test output file
    cortest_output,
    # Imputed data
    imputed_output,
    # peakRAM output
    peakram_output)
  # Stop process from running
  print("There are 0 SNPs in scaffold, stopping...")
}

if (!is.null(ncol(dat$gen))) {
  # Env vars ----------------------------------------------------------------
  
  # Extract and standardize environmental variables and make into dataframe
  env <- terra::extract(dat$envlayers, dat$coords, ID = FALSE)
  env <- scale(env, center = TRUE, scale = TRUE)
  env <- data.frame(env)
  # When only one env layer provided, env colnames will be named simply 'env' which is not informative
  if (ncol(env) == 1) colnames(env) <- names(dat$envlayers)
  
  # Impute any NA values from nearby samples within a given radius
  if (any(is.na(env))) {
    # Count NAs in each column
    have_nas <- colSums(is.na(env))[colSums(is.na(env)) > 0]
    coords_sf <- sf::st_as_sf(dat$coords, coords = c("x", "y"), crs = 4326)

    # Impute to the value of the closest sample
    df <- bind_cols(coords_sf, env)
    mod_sf <- df %>% sf::st_as_sf()

    close_impute <- function(var, mod_sf){
      map_dbl(1:nrow(mod_sf), \(i){ 
        if (!is.na(sf::st_drop_geometry(mod_sf[i, var]))) return(mod_sf[i, var][[1]])
        # Get the closest sample that is not NA
        nearest <- sf::st_nearest_feature(mod_sf[i,], drop_na(mod_sf[-i,]), 1)
        # Get the value of the closest sample
        return(drop_na(mod_sf)[nearest, var][[1]])
      }) 
    }
    vars <- names(have_nas)
    names(vars) <- vars
    sf_imputed <- map(vars, ~close_impute(.x, mod_sf), .progress = TRUE) 
    df_imputed <- bind_cols(dplyr::select(df, -all_of(vars)), bind_cols(sf_imputed))

    # Extract imputed env values for all samples
    env <- df_imputed %>% sf::st_drop_geometry() %>% dplyr::select(colnames(env)) %>% as.data.frame()
  }

  # Run RDA -----------------------------------------------------------------
  
  #' Redefined algatr `rda_run()` function to take in Plink PC files for pRDA correctPC
  #'
  #' @param gen genotype dosage matrix (rows = individuals & columns = SNPs) or `vcfR` object
  #' @param env dataframe with environmental data or a Raster* type object from which environmental values for the coordinates can be extracted
  #' @param coords dataframe with coordinates (only needed if correctGEO = TRUE) or if env is a Raster* from which values should be extracted
  #' @param model whether to fit the model with all variables ("full") or to perform variable selection to determine the best set of variables ("best"); defaults to "full"
  #' @param correctGEO whether to condition on geographic coordinates (defaults to FALSE)
  #' @param correctPC path to eigenvector file produced from Plink --pca to condition on PCs from PCA of genotypes
  #' @param nPC if `correctPC` not NULL, number of PCs to use (defaults to 3)
  #' @param pc_data_type for partial RDA, whether to correct using PCs ("pc") or admixture ("admixture")
  #' @param species_ks if `correctPC = TRUE` and `pc_data_type = "admixture"`, how many Ks to correct
  #' @param Pin if `model = "best"`, limits of permutation P-values for adding (`Pin`) a term to the model, or dropping (`Pout`) from the model. Term is added if` P <= Pin`, and removed if `P > Pout` (see \link[vegan]{ordiR2step}) (defaults to 0.05)
  #' @param R2permutations if `model = "best"`, number of permutations used in the estimation of adjusted R2 for cca using RsquareAdj (see \link[vegan]{ordiR2step}) (defaults to 1000)
  #' @param R2scope if `model = "best"` and set to TRUE (default), use adjusted R2 as the stopping criterion: only models with lower adjusted R2 than scope are accepted (see \link[vegan]{ordiR2step})
  #'
  #' @return RDA model
  #' @export
  rda_run_pc <- function(gen, env, coords = NULL, model = "full", correctGEO = FALSE, correctPC, pc_path, nPC = 3,
                         pc_data_type, species_ks, Pin = 0.05, R2permutations = 1000, R2scope = T) {
    
    # Check that env var names don't match coord names
    if (any(colnames(coords) %in% colnames(env))) {
      colnames(env) <- paste(colnames(env), "_env", sep = "")
      warning("env names should differ from x and y. Appending 'env' to env names")
    }

    if (correctPC) {
      # Read in PCs
      if (pc_data_type == "pca") {
        pc <- readr::read_tsv(paste0(pc_path)) %>%
        tibble::column_to_rownames(var = "#IID") %>% 
        dplyr::select(tidyselect::all_of(1:nPC))
      }
      if (pc_data_type == "admixture") {
        pc <- readr::read_tsv(paste0(pc_path)) %>%
        dplyr::select(-c(likely_assignment, samp_num, sample_order)) %>% 
        tibble::column_to_rownames(var = "Sample") %>% 
        dplyr::select(tidyselect::all_of(1:species_ks))
      }
    }
   
    # Handle NA values -----------------------------------------------------
    if (any(is.na(gen))) {
      stop("Missing values found in gen data")
    }
    
    if (any(is.na(env))) {
      if (correctPC) {
        warning("Missing values found in env data, removing rows with NAs")
        na_env <- env
        gen <- gen[complete.cases(na_env), ]
        pc <- pc[complete.cases(na_env), ]
        # Must come last
        env <- env[complete.cases(na_env), ]
        if (!is.null(coords)) coords <- coords[complete.cases(na_env), ]
        # Check env var naming ----------------------------------------------------
        if(any(colnames(pc) %in% colnames(env))) {
        colnames(env) <- paste(colnames(env), "_env", sep = "")
        warning("env names should differ from PC1, PC2, etc if correctPC is TRUE. Appending 'env' to env names")
    }
      }
      if (!correctPC) {
        warning("Missing values found in env data, removing rows with NAs")
        na_env <- env
        gen <- gen[complete.cases(na_env), ]
        # Must come last
        env <- env[complete.cases(na_env), ]
        if (!is.null(coords)) coords <- coords[complete.cases(na_env), ]
      }
    }
    
    # Set up model ------------------------------------------------------------
    print(paste0("env object has ", nrow(env), " rows"))
    # print(paste0("pc object has ", nrow(pc), " rows"))
    print(paste0("gen object has ", nrow(gen), " rows"))
    
    if (!correctPC) {
      moddf <- data.frame(env)
      f <- as.formula(paste0("gen ~ ", paste(colnames(env), collapse = "+")))
    }
    
    if (correctPC) {
      moddf <- data.frame(env, pc)
      f <- as.formula(paste0("gen ~ ", paste(colnames(env), collapse = "+"), "+ Condition(", paste(colnames(pc), collapse = "+"), ")"))
    }
    
    mod <- vegan::rda(f, data = moddf)
    return(list(gen = gen, mod = mod, env = env))
  }
  
  peakRAM_run <-
    peakRAM::peakRAM(
      mod <- rda_run_pc(gen = dat$gen,
                        env = env,
                        coords = dat$coords, 
                        model = model, 
                        correctGEO = correctGEO, 
                        correctPC = correctPC, 
                        pc_path = pc_path,
                        nPC = nPC, 
                        pc_data_type = pc_data_type,
                        species_ks = species_ks,
                        Pin = Pin, 
                        R2permutations = R2permutations, 
                        R2scope = R2scope)
    )
  
  
  # Get outliers and run cortest --------------------------------------------
  
  # Determine RDA outliers based on Z-scores; function from algatr package
  z_outlier_method <- function(mod, naxes, z) {
    if (naxes == "all") naxes <- ncol(mod$CCA$v)
    load.rda <- vegan::scores(mod, choices = 1:naxes, display = "species")
    
    results <-
      purrr::map(1:ncol(load.rda), ~z_outlier_helper(.x, load.rda, z)) %>%
      dplyr::bind_rows()
    
    return(results)
  }
  # z_outlier_method helper function; function from algatr package
  z_outlier_helper <- function(axis, load.rda, z) {
    x <- load.rda[, axis]
    out <- outliers(x, z)
    cand <- cbind.data.frame(names(out), rep(axis, times = length(out)), unname(out))
    colnames(cand) <- c("rda_snps", "axis", "loading")
    cand$rda_snps <- as.character(cand$rda_snps)
    return(cand)
  }
  # Z-scores outlier finder
  outliers <- function(x, z) {
    lims <- mean(x) + c(-1, 1) * z * sd(x) # find loadings +/-z sd from mean loading
    x[x < lims[1] | x > lims[2]] # SNP names in these tails
  }
  
  get_outliers <- function(mod, gen) {
    safe_rda_getoutliers <- purrr::safely(algatr::rda_getoutliers)
    safe_z_outlier_method <- purrr::safely(z_outlier_method)
    # safe_rda_sig_z <- safe_rda_getoutliers(mod, naxes = "all", outlier_method = "z", z = z, plot = FALSE)
    safe_rda_sig_z <- safe_z_outlier_method(mod, naxes = "all", z = z)
    #safe_rda_sig_p <- safe_rda_getoutliers(mod, naxes = "all", outlier_method = "p", p_adj = p_adj, sig = sig, plot = FALSE)
    
    # If any of the above error out, assign NULL value to object in result
    if (is.null(safe_rda_sig_z$error)) {
      if (nrow(safe_rda_sig_z$result) > 0) {
        rda_sig_z <- safe_rda_sig_z$result
        # Extract genotypes for outlier SNPs
        rda_snps_z <- rda_sig_z$rda_snps
        rda_gen_z <- gen[, rda_snps_z]
      }
      # If no significant SNPs found:
      if (nrow(safe_rda_sig_z$result) == 0) {
        print("No significant SNPs found using Z-scores method")
        rda_sig_z <- NULL
        rda_gen_z <- NULL
      }
    }
    if (!is.null(safe_rda_sig_z$error)) {
      rda_sig_z <- NULL
      rda_gen_z <- NULL
      print(safe_rda_sig_z$error)
    }
    
#    if (is.null(safe_rda_sig_p$error)) {
 #     if (length(safe_rda_sig_p$result$rda_snps) > 0) {
  #      rda_sig_p <- safe_rda_sig_p$result
   #     # Extract genotypes for outlier SNPs
    #    rda_snps_p <- rda_sig_p$rda_snps
     #   rda_gen_p <- gen[, rda_snps_p]
      #}
      ## If no significant SNPs found:
      #if (length(safe_rda_sig_p$result$rda_snps) == 0) {
       # print("No significant SNPs found using p-values method")
        #rda_sig_p <- NULL
        #rda_gen_p <- NULL
      #}
      ## p-values can't be calculated if there are fewer than 2 RDA axes
      #if (ncol(mod$CCA$v) < 2) {
       # print("Skipping p-value outlier detection; fewer than 2 RDA axes")
        #rda_sig_p <- NULL
        #rda_gen_p <- NULL 
 #     }
  #  }
   # if (!is.null(safe_rda_sig_p$error)) {
    #  rda_sig_p <- NULL
     # rda_gen_p <- NULL
      #print(safe_rda_sig_p$error)
    #}
    
    rda_sig_p <- NULL
    rda_gen_p <- NULL

    results <- list(rda_sig_z = rda_sig_z,
                    rda_sig_p = rda_sig_p,
                    rda_gen_p = rda_gen_p,
                    rda_gen_z = rda_gen_z)
    
    return(results)
  }
  
  peakRAM_outliers <-
    peakRAM::peakRAM(
      outliers <- get_outliers(mod = mod$mod, gen = mod$gen)
    )
  
  print(summary(outliers)) # to assess whether there are NULL objects
  
  # Run correlation test
  run_cortest <- function(rda_gen_p, rda_gen_z, env) {
    # Check if either outlier detection method is NULL
    if (is.null(rda_gen_p)) cor_df_p <- NULL
    if (is.null(rda_gen_z)) cor_df_z <- NULL
    
    # If outliers detected, run correlation test
    if (!is.null(rda_gen_p)) cor_df_p <- algatr::rda_cor(rda_gen_p, env)
    if (!is.null(rda_gen_z)) cor_df_z <- algatr::rda_cor(rda_gen_z, env)
    
    results <- list(cor_df_p = cor_df_p, 
                    cor_df_z = cor_df_z)
    
    return(results)
  }
  
  peakRAM_cortest <-
    peakRAM::peakRAM(
      cortest <- run_cortest(rda_gen_p = outliers$rda_gen_p, rda_gen_z = outliers$rda_gen_z, env = mod$env)
    )
  print(summary(cortest))
  
  # Export results ----------------------------------------------------------
  
  export_rda <- function(mod, rda_sig_z, rda_sig_p, cor_df_p, cor_df_z, save_impute) {
    outlier_helper <- function(df, outlier) {dat <- df %>% dplyr::mutate(outlier_method = outlier)}
    
    # RDA model results -------------------------------------------------------
    # Export raw values from RDA model result
    data.frame(mod$colsum) %>% rownames_to_column(var = "locus") %>% write_csv(file = colsum_output, col_names = TRUE)
    data.frame(mod$Ybar) %>% rownames_to_column(var = "INDV") %>% write_csv(file = ybar_output, col_names = TRUE)
    data.frame(mod$CCA$v) %>% rownames_to_column(var = "locus") %>% write_csv(file = v_output, col_names = TRUE)
    data.frame(mod$CCA$u) %>% rownames_to_column(var = "INDV") %>% write_csv(file = u_output, col_names = TRUE)
    data.frame(mod$CCA$wa) %>% rownames_to_column(var = "INDV") %>% write_csv(file = wa_output, col_names = TRUE)
    data.frame(mod$CCA$QR$qr) %>% write_csv(file = qr_output, col_names = TRUE)
    data.frame(mod$CCA$eig) %>% rownames_to_column(var = "RDA") %>% write_csv(file = eig_output, col_names = TRUE)
    data.frame(mod$CCA$biplot) %>% rownames_to_column(var = "var") %>% write_csv(file = biplot_output, col_names = TRUE)
    data.frame(mod$CCA$QR$qraux) %>% write_csv(file = qraux_output, col_names = TRUE)
    data.frame(mod$CCA$envcentre) %>% tibble::rownames_to_column(var = "axis") %>% write_csv(file = envcentre_output, col_names = TRUE)
    data.frame(mod_chi = mod$tot.chi,
               mod_chi_cca = mod$CCA$tot.chi) %>% write_csv(file = chi_output, col_names = TRUE)
    # Export scores (scaled and unscaled); labeled "species" corresponds to v, "sites" corresponds to wa, and "constraints" corresponds to u
    scaled_loadings <- vegan::scores(mod, choices = 1:ncol(mod$CCA$v), tidy = TRUE) %>% write_csv(file = scalload_output, col_names = TRUE)
    unscaled_loadings <- vegan::scores(mod, choices = 1:ncol(mod$CCA$v), tidy = TRUE, scaling = 0) %>% write_csv(file = unscalload_output, col_names = TRUE)

    # Outlier results ---------------------------------------------------------
    # Sig results Z-scores
    if (!is.null(rda_sig_z)) readr::write_csv(rda_sig_z,
                                              #file = paste0(output_path, species, "_RDA_outliers_", model, "_Zscores.csv"),
                                              file = zscore_output,
                                              col_names = TRUE)
    # Sig results p-values
    if (!is.null(rda_sig_p)) {
      snps <- rda_sig_p$rdadapt %>% 
        dplyr::mutate(locus = colnames(dat$gen))
      readr::write_csv(snps,
                       file = rdadapt_output,
                       col_names = TRUE)
    }
    
    # If NULL results for outliers
    if (is.null(rda_sig_z)) file.create(zscore_output)
    if (is.null(rda_sig_p)) file.create(rdadapt_output)
    
    # Correlation test results ------------------------------------------------
    if (!is.null(cor_df_p) | !is.null(cor_df_z)) {
      if (!is.null(cor_df_p) & !is.null(cor_df_z)) {
        cor_test <- rbind(outlier_helper(cor_df_p, outlier = "p"),
                          outlier_helper(cor_df_z, outlier = "z"))
      }
      if (!is.null(cor_df_p) & is.null(cor_df_z)) cor_test <- outlier_helper(cor_df_p, outlier = "p")
      if (!is.null(cor_df_z) & is.null(cor_df_p)) cor_test <- outlier_helper(cor_df_z, outlier = "z")
      readr::write_csv(cor_test, 
                       file = cortest_output,
                       col_names = TRUE)
    }
    
    if (is.null(cor_df_z) & is.null(cor_df_p)) file.create(cortest_output)
    
    # Save imputed data -------------------------------------------------------
    if (save_impute) {
      write.table(dat$gen, 
                  file = imputed_output,
                  sep = " ", row.names = FALSE, col.names = TRUE, quote = FALSE)
    }
  
  # Export results
  peakRAM_exp <- 
    peakRAM::peakRAM(
      export_rda(mod$mod, outliers$rda_sig_z, outliers$rda_sig_p, cortest$cor_df_p, cortest$cor_df_z, save_impute = save_impute)
    )
  
  # Export RAM usage --------------------------------------------------------
  
  RAM <- dplyr::bind_rows(peakRAM_imp,
                          peakRAM_run, 
                          peakRAM_outliers,
                          peakRAM_cortest,
                          peakRAM_exp) %>% 
    dplyr::mutate(fxn = c("import", "run", "outliers", "cortest", "export"))
  
  readr::write_csv(RAM,
                   file = peakram_output)
}

