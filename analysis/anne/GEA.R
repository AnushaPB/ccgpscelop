#' Redefined algatr `rda_run()` function to take in Plink PC files for pRDA correctPC
#'
#' @param gen genotype dosage matrix (rows = individuals & columns = SNPs) or `vcfR` object
#' @param env dataframe with environmental data or a Raster* type object from which environmental values for the coordinates can be extracted
#' @param coords dataframe with coordinates (only needed if correctGEO = TRUE) or if env is a Raster* from which values should be extracted
#' @param model whether to fit the model with all variables ("full") or to perform variable selection to determine the best set of variables ("best"); defaults to "full"
#' @param correctGEO whether to condition on geographic coordinates (defaults to FALSE)
#' @param correctPC path to eigenvector file produced from Plink --pca to condition on PCs from PCA of genotypes; defaults to NULL
#' @param nPC if `correctPC` not NULL, number of PCs to use (defaults to 3)
#' @param Pin if `model = "best"`, limits of permutation P-values for adding (`Pin`) a term to the model, or dropping (`Pout`) from the model. Term is added if` P <= Pin`, and removed if `P > Pout` (see \link[vegan]{ordiR2step}) (defaults to 0.05)
#' @param R2permutations if `model = "best"`, number of permutations used in the estimation of adjusted R2 for cca using RsquareAdj (see \link[vegan]{ordiR2step}) (defaults to 1000)
#' @param R2scope if `model = "best"` and set to TRUE (default), use adjusted R2 as the stopping criterion: only models with lower adjusted R2 than scope are accepted (see \link[vegan]{ordiR2step})
#'
#' @return RDA model
#' @export

rda_run_pc <- function(gen, env, coords = NULL, model = "full", correctGEO = FALSE, correctPC = NULL, nPC = 3,
                       Pin = 0.05, R2permutations = 1000, R2scope = T) {

  # Check that env var names don't match coord names
  if (any(colnames(coords) %in% colnames(env))) {
    colnames(env) <- paste(colnames(env), "_env", sep = "")
    warning("env names should differ from x and y. Appending 'env' to env names")
  }

  # Read in PCs
  pc <- readr::read_tsv(paste0(correctPC)) %>%
    tibble::column_to_rownames(var = "#IID") %>%
    dplyr::select(tidyselect::all_of(1:nPC))

  # Handle NA values -----------------------------------------------------
  if (any(is.na(gen))) {
    stop("Missing values found in gen data")
  }

  if (any(is.na(env))) {
    warning("Missing values found in env data, removing rows with NAs")
    na_env <- env
    gen <- gen[complete.cases(na_env), ]
    pc <- pc[complete.cases(na_env), ]
    # Must come last
    env <- env[complete.cases(na_env), ]
    if (!is.null(coords)) coords <- coords[complete.cases(na_env), ]
  }

  # Set up model ------------------------------------------------------------
  # Check env var naming ----------------------------------------------------
  if(any(colnames(pc) %in% colnames(env))) {
    colnames(env) <- paste(colnames(env), "_env", sep = "")
    warning("env names should differ from PC1, PC2, etc if correctPC is TRUE. Appending 'env' to env names")
  }
  print(paste0("env object has ", nrow(env), " rows"))
  print(paste0("pc object has ", nrow(pc), " rows"))
  print(paste0("gen object has ", nrow(gen), " rows"))

  moddf <- data.frame(env, pc)
  f <- as.formula(paste0("gen ~ ", paste(colnames(env), collapse = "+"), "+ Condition(", paste(colnames(pc), collapse = "+"), ")"))

  mod <- vegan::rda(f, data = moddf)
  return(list(gen = gen, mod = mod, env = env))
}

#' Get outliers using both p-value and Z-score methods
#'
#' @param mod RDA model
#' @param gen genetic data
#' @param p_adj
#' @param sig
#'
#' @return
#' @export
get_outliers <- function(mod, gen, p_adj, sig) {
  rda_sig_z <- rda_getoutliers(mod, naxes = "all", outlier_method = "z", z = 3, plot = FALSE)
  rda_sig_p <- rda_getoutliers(mod, naxes = "all", outlier_method = "p", p_adj = p_adj, sig = sig, plot = FALSE)

  # Extract genotypes for outlier SNPs
  rda_snps_p <- rda_sig_p$rda_snps
  rda_gen_p <- gen[, rda_snps_p]
  rda_snps_z <- rda_sig_z$rda_snps
  rda_gen_z <- gen[, rda_snps_z]

  results <- list(rda_sig_z = rda_sig_z,
                  rda_sig_p = rda_sig_p,
                  rda_gen_p = rda_gen_p,
                  rda_gen_z = rda_gen_z)

  return(results)
}

#' Run RDA correlation tests
#'
#' @param rda_gen_p
#' @param rda_gen_z
#' @param env
#'
#' @return
#' @export
run_cortest <- function(rda_gen_p, rda_gen_z, env) {
  # Run correlation test
  cor_df_p <- algatr::rda_cor(rda_gen_p, env)
  cor_df_z <- algatr::rda_cor(rda_gen_z, env)

  results <- list(cor_df_p = cor_df_p,
                  cor_df_z = cor_df_z)

  return(results)
}

#' Make Manhattan plot of RDA results
#'
#' @param mod
#' @param outliers
#'
#' @return
#' @export
manhat_plot <- function(mod, outliers) {
  # Make and get tidy data frames for plotting
  snp_scores <- vegan::scores(mod, choices = 1:ncol(mod$CCA$v), display = "species", scaling = "none")
  TAB_snps <- data.frame(names = row.names(snp_scores), snp_scores)
  TAB_snps$type <- "Non-outlier"
  TAB_snps$type[TAB_snps$names %in% outliers$rda_sig_p$rda_snps] <- "Outlier"
  TAB_snps$type <- factor(TAB_snps$type, levels = c("Non-outlier", "Outlier"))

  TAB_manhattan <- data.frame(
    pos = 1:nrow(TAB_snps),
    pvalues = outliers$rda_sig_p$pvalues,
    type = factor(TAB_snps$type, levels = c("Non-outlier", "Outlier"))
  )

  TAB_manhattan <- TAB_manhattan %>%
    tibble::rownames_to_column(var = "name") %>%
    tidyr::separate_wider_delim(cols = name,
                                delim = "_",
                                names = c("chrom", "site"))

  TAB_manhattan <- TAB_manhattan[order(TAB_manhattan$pos), ]
  TAB_manhattan$chrom <- factor(TAB_manhattan$chrom, levels = (unique(TAB_manhattan$chrom)))

  ylim <- TAB_manhattan %>%
    dplyr::filter(pvalues == min(pvalues)) %>%
    mutate(ylim = abs(floor(log10(pvalues))) + 2) %>%
    pull(ylim)

  axis_set <- TAB_manhattan %>%
    group_by(chrom) %>%
    summarize(center = mean(pos))

  plt_manhat <-
    ggplot2::ggplot() +
    ggplot2::geom_point(data = TAB_manhattan %>% dplyr::filter(type == "Outlier"),
                        ggplot2::aes(x = pos, y = -log10(pvalues),), col = "orange", size = 1.4, alpha = 0.75) +
    ggplot2::xlab(NULL) +
    ggplot2::ylab("-log10(p)") +
    ggplot2::geom_hline(yintercept = -log10(sig), linetype = "dashed", color = "black", linewidth = 0.6) +
    scale_x_continuous(label = axis_set$chrom, breaks = axis_set$center) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, ylim)) +
    ggplot2::geom_point(data = TAB_manhattan %>% dplyr::filter(type == "Non-outlier"),
                        ggplot2::aes(x = pos, y = -log10(pvalues), col = chrom), size = 1.4, alpha = 0.75) +
    ggplot2::scale_color_manual(values = rep(c("#276FBF","#183059"),
                                             ceiling(length(unique(TAB_manhattan$chrom))/2))[1:length(unique(TAB_manhattan$chrom))]) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      legend.position = "none",
      panel.grid = ggplot2::element_blank(),
      plot.background = ggplot2::element_blank(),
      axis.text.x = element_text(angle = 60, size = 4, vjust = 0.5)
    )

  return(plt_manhat)
}

#' Export RDA results
#'
#' @param mod
#' @param output_path
#' @param species
#' @param model
#' @param rda_sig_z
#' @param rda_sig_p
#' @param cor_df_p
#' @param cor_df_z
#' @param save_impute
#'
#' @return
#' @export
export_rda <- function(mod, output_path, species, model, rda_sig_z, rda_sig_p, cor_df_p, cor_df_z, save_impute) {
  outlier_helper <- function(df, outlier) {dat <- df %>% dplyr::mutate(outlier_method = outlier)}

  # RDA model results
  saveRDS(mod, file = paste0(output_path, species, "_RDA_model_", model, ".RDS"))

  # Sig results Z-scores
  readr::write_csv(rda_sig_z,
                   file = paste0(output_path, species, "_RDA_outliers_", model, "_Zscores.csv"),
                   col_names = TRUE)

  # Sig results p-values
  snps <- rda_sig_p$rdadapt %>%
    dplyr::mutate(locus = colnames(dat$gen))
  readr::write_csv(snps,
                   file = paste0(output_path, species, "_RDA_outliers_", model, "_rdadapt.csv"),
                   col_names = TRUE)

  # Correlation test results
  cor_test <- rbind(outlier_helper(cor_df_p, outlier = "p"),
                    outlier_helper(cor_df_z, outlier = "z"))
  readr::write_csv(cor_test, file = paste0(output_path, species, "_RDA_cortest_", model, ".csv"),
                   col_names = TRUE)

  # Save imputed data
  if (save_impute) {
    write.table(dat$gen, file = paste0(output_path, species, "_imputed_", impute, ".txt"),
                sep = " ", row.names = FALSE, col.names = TRUE, quote = FALSE)
  }

  # Build and save Manhattan plot to file
  #plt_manhat <- manhat_plot(mod, outliers)
  #plt_manhat
  #ggsave(paste0(output_path, species, "_manhatplot.png"), width = 8, height = 4.5, bg = "white")
}

