library(here)
library(tidyverse)
library(vcfR)
library(algatr)
library(raster)
library(tidyterra)
library(terra)
library(viridis)
library(cowplot)

### Plotting genotypes of topmost outliers across landscape
# Run within analysis/anne directory

# Read in files -----------------------------------------------------------

rdadapt <- read_csv(here("outputs/RDA/58-Sceloporus_RDA_outliers_full_rdadapt.csv"))

# Look at distribution of SNPs
rdadapt %>% 
    # filter(q.values <= 0.01) %>% 
    ggplot(aes(x = q.values)) +
    geom_histogram()
ggsave(here("analysis/anne/outputs/RDA_allSNPs.pdf")) # allsig is filtered at 0.01

# Pull out topmost outliers for chr1
sig = 0.0000000000000000000000000001 # 1x10e-28
chr = "chr1"

# 296 SNPs
supersig <-
    rdadapt %>% 
    filter(scaff == chr) %>% 
    filter(q.values <= sig)
    
write_tsv(supersig %>% 
    tidyr::extract(locus, c("CHROM", "POS", "REF", "ALT"), regex = "(.*)_([^_]+)_([^_]+)_([^_]+)$") %>% 
    dplyr::select(CHROM, POS), 
    here("analysis/anne/outputs/supersig_SNPs.txt"), col_names = FALSE)

# Extract relevant SNPs from vcf in bash using genotype_dists_gea.sh script

# Read in supersig vcf to retrieve genotypes
supersig_vcf <- read.vcfR(here("analysis/anne/outputs/58-Sceloporus_chr1_supersigSNPs.vcf"))
dos <- vcf_to_dosage(supersig_vcf)

env <- terra::rast(here("data/env/california_chelsa_bioclim_1981-2010_V.2.1_pca.tif"))
coords <- read_tsv(here("data/ccgp_data/58-Sceloporus.coords.txt"), col_names = c("INDV", "x", "y"))

genID <- rownames(dos)
overlap <- coords$INDV %in% genID
coordsF <- coords[overlap,]
# Check order
coordsF <- coordsF[match(genID, coordsF$INDV),]

supersig <- supersig %>% filter(locus %in% colnames(dos))


# Build map -----------------------------------------------------------

# Top 5 SNPs
topsnps <- supersig %>% slice_max(q.values, n = 5)

subset <- as.data.frame(dos) %>% dplyr::select(topsnps$locus) %>% rownames_to_column(var = "INDV") %>% pivot_longer(cols = 2:6, names_to = "locus", values_to = "GT")
subset <- left_join(subset, coordsF)
subset$GT <- as.factor(subset$GT)

# Facet wrap on locus
ggplot() +
    geom_spatraster(data = env[[1]]) +
    scale_fill_viridis_c(option = "D", na.value = "transparent") +
    theme_map() +
    geom_point(data = subset %>% na.omit(), aes(x = x, y = y, color = GT)) +
    scale_color_manual(values = c("yellow", "orange", "red")) +
    ggtitle(names(env[[1]])) +
    facet_wrap(~locus)
ggsave(here("analysis/anne/outputs/supersigsnps_envPC1.pdf"))
ggplot() +
    geom_spatraster(data = env[[2]]) +
    scale_fill_viridis_c(option = "D", na.value = "transparent") +
    theme_map() +
    geom_point(data = subset %>% na.omit(), aes(x = x, y = y, color = GT)) +
    scale_color_manual(values = c("yellow", "orange", "red")) +
    ggtitle(names(env[[2]])) +
    facet_wrap(~locus)
ggsave(here("analysis/anne/outputs/supersigsnps_envPC2.pdf"))
ggplot() +
    geom_spatraster(data = env[[3]]) +
    scale_fill_viridis_c(option = "D", na.value = "transparent") +
    theme_map() +
    geom_point(data = subset %>% na.omit(), aes(x = x, y = y, color = GT)) +
    scale_color_manual(values = c("yellow", "orange", "red")) +
    ggtitle(names(env[[3]])) +
    facet_wrap(~locus)
ggsave(here("analysis/anne/outputs/supersigsnps_envPC3.pdf"))
# plot_grid(pc1, pc2, pc3, ncol = 3)
# ggsave(here("analysis/anne/outputs/envPCs_supersigsnps.pdf"))

# Facet on env var
plot_genotypes <- function(envlayer, subset, locus) {
    ggplot() +
        geom_spatraster(data = envlayer) +
        scale_fill_viridis_c(option = "D", na.value = "transparent") +
        theme_map() +
        geom_point(data = subset %>% filter(locus == locus) %>% na.omit(), aes(x = x, y = y, color = GT)) +
        scale_color_manual(values = c("yellow", "orange", "red")) +
        ggtitle(paste0(locus, ", env ", names(envlayer)))
}
pc1 <- plot_genotypes(envlayer = env[[1]], subset, locus = chr)
pc2 <- plot_genotypes(envlayer = env[[2]], subset, locus = chr)
pc3 <- plot_genotypes(envlayer = env[[3]], subset, locus = chr)
plot_grid(pc1, pc2, pc3, ncol = 3)
ggsave(here("analysis/anne/outputs/testpc.pdf"))

# # Try to plot a different way
# pt_dat <- terra::as.data.frame(env[[1]], xy = TRUE, na.rm = TRUE)
# pt_dat$x <- round(pt_dat$x, 2)

# pt_dat <- left_join(pt_dat, subset)

# ggplot(pt_dat %>% na.omit()) +
#     geom_point(aes(x = GT, y = PC1))
#     # geom_histogram(aes(x = PC1))
# ggsave(here("analysis/anne/outputs/test.pdf"))
