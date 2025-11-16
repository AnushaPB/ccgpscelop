library(tidyverse)
library(here)
library(sf)
source(here("general_functions.R"))
dos <- read_table(here("analysis", "check_nonsyn", "outputs", "all_nonsynonymous_dosage.raw"), col_names = TRUE)
coords <- get_coords(sf = TRUE)
dos <- dos %>% rename(SampleID = IID) %>% filter(SampleID %in% coords$SampleID)
coords <- coords %>% filter(SampleID %in% dos$SampleID)

library(sf)
library(dplyr)
library(purrr)

# coords is already sf = TRUE and in a projected CRS (NAD83 / California Albers)
# dos is your tibble with first 6 PLINK columns + SNP columns

# 1. Distance matrix between all individuals ---------------------------------

# (Assume rows of `coords` and `dos` are in the same order!)
stopifnot(nrow(coords) == nrow(dos))

# pairwise distances (returned in CRS units, here meters)
D <- st_distance(coords)                  # n x n units matrix
D_mat <- as.matrix(D)                     # drop units for speed, keep meters

# 2. Pull genotype matrix (just SNP columns) ---------------------------------
geno_mat <- as.matrix(dos[ , -(1:6)])     # drop FID, IID, PAT, MAT, SEX, PHENOTYPE
n_ind  <- nrow(geno_mat)
n_snp  <- ncol(geno_mat)

# 3. Function to compute "other-allele distance" for one SNP -----------------

other_allele_dist_one_snp <- function(g, D_mat) {
  # g: numeric vector of length n_ind with 0,1,2,NA
  out <- rep(NA_real_, length(g))

  missing <- is.na(g)
  out[missing] <- NA_real_

  # heterozygotes: have both alleles locally
  het <- which(g == 1)
  out[het] <- 0

  # homozygous reference: look for nearest 1 or 2 (alt present)
  ref  <- which(g == 0)
  alt_carriers <- which(g %in% c(1, 2))

  if (length(ref) > 0 && length(alt_carriers) > 0) {
    out[ref] <- apply(D_mat[ref, alt_carriers, drop = FALSE], 1, min)
  }

  # homozygous alt: look for nearest 0 or 1 (ref present)
  alt_hom <- which(g == 2)
  ref_carriers <- which(g %in% c(0, 1))

  if (length(alt_hom) > 0 && length(ref_carriers) > 0) {
    out[alt_hom] <- apply(D_mat[alt_hom, ref_carriers, drop = FALSE], 1, min)
  }

  out
}

# 4. Apply across all SNPs ----------------------------------------------------

# This will return an n_ind x n_snp matrix of distances (in meters)
library(furrr)
plan(multisession, workers = 4)
dist_mat <- future_map_dfc(
  .x = as.data.frame(geno_mat),
  .f = ~ other_allele_dist_one_snp(.x, D_mat),
  .progress = TRUE
)
plan(sequential)

# Take average across all columns (SNPs)
dist_avg <- rowMeans(dist_mat, na.rm = TRUE)
coords$dist_to_allele <- dist_avg


pdf(here("analysis", "check_nonsyn", "plots", "avg_dist_to_allele.pdf"), width = 6, height = 4)
ggplot(coords) +
  geom_sf(aes(color = dist_to_allele), size = 2) +
  scale_color_viridis_c(option = "plasma", na.value = "grey80") +
  labs(color = "Avg. distance to\nother allele (m)") +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_line(color = "transparent"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )
dev.off()
