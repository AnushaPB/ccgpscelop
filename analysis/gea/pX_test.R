library(tidyverse)
library(here)
dos <- read_table(here("analysis", "gea", "outputs", "thinned.raw"))
gea_dos <- read_table(here("analysis", "gea", "outputs", "genes.raw"))
pop_df <- read_csv(here("analysis", "admixture", "outputs", "Q9.csv")) # socal: cluster = 9

# Load coords and model dataframe
coords <- 
  get_coords(sf = TRUE) %>% 
  filter(SampleID %in% dos$IID) %>% 
  st_transform(3310) %>%
  # REMOVE CHANNEL ISLAND
  filter(SampleID != "Scelocci_CCGPMC_MW01-3-14") %>%
  left_join(pop_df) %>%
  mutate(socal = as.numeric(cluster == 9))

mod_df <- read_csv(here("analysis", "genetic_diversity", "outputs", "model_df.csv"))
mod_df <- left_join(coords, mod_df, by = "SampleID")

# Subset to only samples in coords
dos <- dos[dos$IID %in% coords$SampleID, ]
# Confirm orders are the sam
stopifnot(all(dos$IID == coords$SampleID))
# Convert to heterozygosity using ifelse for matrices
dosmat <- as.matrix(dos[, -(1:6)])
het <- matrix(ifelse(dosmat == 1, 1, 0), nrow = nrow(dosmat), ncol = ncol(dosmat))
colnames(het) <- colnames(dosmat)

gea <- gea_dos[gea_dos$IID %in% coords$SampleID, ]
# Confirm orders are the same
stopifnot(all(gea$IID == coords$SampleID))
# Convert to heterozygosity using ifelse for matrices
gea_mat <- as.matrix(gea[, -(1:6)])
het_gea <- matrix(ifelse(gea_mat == 1, 1, 0), nrow = nrow(gea_mat), ncol = ncol(gea_mat))

# for each column of dos calculate the correlation with 
corp <- function(x, y){
  cor.test(x, y, method = "pearson", use = "pairwise.complete.obs")$p.value
}

# Apply corp across cols of dosmat
dos_socal <- map_dbl(1:ncol(dosmat), ~corp(x = coords$socal, y = dosmat[, .x]), .progress = TRUE)
names(dos_socal) <- colnames(dosmat)
socal_sig <- (dos_socal < 0.1)
mod_df$socal_het <- rowMeans(het[, socal_sig], na.rm = TRUE)
mod_df$notsocal_het <- rowMeans(het[, !socal_sig], na.rm = TRUE)

dos_bio1 <- map_dbl(1:ncol(dosmat), ~corp(x = mod_df$bio1, y = dosmat[, .x]), .progress = TRUE)
names(dos_bio1) <- colnames(dosmat)
bio1_sig <- (dos_bio1 < 0.05)
mod_df$het_bio1 <- rowMeans(het[, bio1_sig], na.rm = TRUE)

gea_socal <- map_dbl(1:ncol(gea_mat), ~corp(x = coords$socal, y = gea_mat[, .x]), .progress = TRUE)
names(gea_socal) <- colnames(gea_mat)
socal_sig_gea <- (gea_socal < 0.1)
mod_df$socal_het_gea <- rowMeans(het_gea[, socal_sig_gea], na.rm = TRUE)
mod_df$notsocal_het_gea <- rowMeans(het_gea[, !socal_sig_gea], na.rm = TRUE)
mod_df$test_gea <- rowMeans(het_gea, na.rm = TRUE)

mean(socal_sig, na.rm = TRUE)
mean(socal_sig_gea, na.rm = TRUE)

create_plot <- function(ycol, xcol = "Ho") {
  ggplot(mod_df, aes(x = .data[[xcol]], y = .data[[ycol]], col = factor(socal))) +
    geom_point() +
    geom_smooth(method = "lm") +
    theme_classic() 
}

plots <- list(
  create_plot("test_gea"),
  create_plot("socal_het_gea"),
  create_plot("notsocal_het_gea"),
  create_plot("socal_het"),
  create_plot("het_bio1"),
  create_plot("notsocal_het")
)


pdf(here("analysis", "gea", "plots", "socal_het.pdf"), width = 20, height = 10)
cowplot::plot_grid(plotlist = plots)
dev.off()

# RDA test
library(algatr)
dosmat_imputed <- simple_impute(dosmat)

mod_full <- rda_run(dosmat_imputed[,1:500000], st_drop_geometry(mod_df[,c("bio1", "NDVI")]), model = "full")
rda_sig_p <- rda_getoutliers(mod_full, naxes = "all", outlier_method = "p", p_adj = "fdr", sig = 0.01, plot = FALSE)
mod_df$rda_het <- rowMeans(het[,1:500000][,rda_sig_p$rda_snps], na.rm = TRUE)

pdf(here(plotpath, "socal_het.pdf"))
mod_pca <- rda_run(dosmat_imputed[,1:500000], st_drop_geometry(mod_df[,c("bio1", "NDVI")]), model = "full", correctPC = TRUE, nPC = 3)
dev.off()
rda_sig_pca <- rda_getoutliers(mod_pca, naxes = "all", outlier_method = "p", p_adj = "fdr", sig = 0.01, plot = FALSE)
mod_df$rdapca_het <- rowMeans(het[,1:500000][,rda_sig_pca$rda_snps], na.rm = TRUE)

pca <- prcomp(dosmat_imputed[,1:100000])
mod_df$PC1 <- pca$x[,1]

mod_pca <- 
plots <- list(
  create_plot("test_gea"),
  create_plot("socal_het_gea"),
  create_plot("notsocal_het_gea"),
  create_plot("socal_het"),
  create_plot("het_bio1"),
  create_plot("notsocal_het"),
  #create_plot("rda_het"),
  #create_plot("rdapca_het"),
  create_plot("rda_het") + scale_y_continuous(limits = c(NA, 0.3)),
  create_plot("rdapca_het") + scale_y_continuous(limits = c(NA, 0.3)),
  create_plot("PC1", xcol = "bio1") + scale_y_continuous(limits = c(NA, 0.3))
)


pdf(here("analysis", "gea", "plots", "socal_het.pdf"), width = 12, height = 10)
cowplot::plot_grid(plotlist = plots)
dev.off()
