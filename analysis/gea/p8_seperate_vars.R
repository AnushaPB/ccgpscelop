library(tidyverse)
library(here)
outpath <- here("analysis", "gea", "outputs")
genes <- read_csv(here(outpath, "bio1ndvi_gea_genes_snp.csv"))

chrs <- paste0("chr", 1:10)
names(chrs) <- chrs

cortest <- 
  map(chrs, ~ read_csv(here(outpath, "RDA_results", .x, "58-Sceloporus_RDA_cortest_full.csv"))) %>%
  bind_rows(.id = "chr") %>%
  dplyr::rename(locus = snp) %>%
  filter(locus %in% genes$locus) %>%
  filter(outlier_method == "z") %>%
  dplyr::select(r, p, var, locus)

# NOTE FOR Z-SCORES ROWS ARE DUPLICATED (MAYBE BECAUSE TWO RDA AXES? ASK ANNE)
cortest <- distinct(cortest)

# Note: some loci in genes will not be in cortest because they are the linked loci
# Conversely, some loci in cortest will not be in genes because they are not in genes
rda_r <- 
  cortest %>%
  # Need to filter before grouping
  filter(p < 0.05) %>%
  dplyr::select(r, var, locus) %>%
  pivot_wider(names_from = var, values_from = r) %>%
  mutate(group = case_when(
    !is.na(`CHELSA_bio1_1981.2010_V.2.1`) & is.na(NDVI) ~ "bio1",
    is.na(`CHELSA_bio1_1981.2010_V.2.1`) & !is.na(NDVI) ~ "ndvi",
    !is.na(`CHELSA_bio1_1981.2010_V.2.1`) & !is.na(NDVI) ~ "both"
  ))

rda_r %>% group_by(group) %>% count()

# Bring in the info about linked loci
linked_info <- 
  read_csv(here(outpath, "bio1ndvi_rda_linked_snps_info.csv")) %>%
  dplyr::select(outlier_locus, linked_locus) %>%
  # Note: some outlier loci do not have linked loci
  drop_na(linked_locus)

joined <- 
  left_join(linked_info, rda_r, by = c("outlier_locus" = "locus")) %>%
  # Filter to loci in genes
  dplyr::filter(linked_locus %in% genes$locus & outlier_locus %in% genes$locus) %>%
  # Drop outlier_locus column and rename linked_locus to locus and bind to rda_r
  dplyr::select(-outlier_locus) %>%
  dplyr::rename(locus = linked_locus) %>%
  bind_rows(rda_r) 

write_csv(joined, here(outpath, "bio1ndvi_genes_cortest.csv"))

# FIGURING OUT MYSTERY (MIGHT BE RESOLVED COME BACK TO THIS LATER)
library(here)
library(tidyverse)
r <- read_table(here("analysis", "gea", "outputs", "snp_r2.ld"))
rdasig <- read_csv(here(outpath, "RDA_bio1_ndvi", "58-Sceloporus_RDA_outliers_full_rdadapt.csv")) 

# Filter pairs where both SNP_A and SNP_B are in rdasig$locus
# (e..g, pairs with R2 > 0.6 where both SNPs ended up in RDA)
filtered_r <- 
  r %>%
  filter(SNP_A %in% rdasig$locus & SNP_B %in% rdasig$locus)

# Count number of pairs
nrow(filtered_r)

# Count number of SNPs
length(unique(c(filtered_r$SNP_A, filtered_r$SNP_B)))

# Example:
snpa <- filtered_r[1, ] %>% pull(SNP_A)
snpb <- filtered_r[1, ] %>% pull(SNP_B)
rdasig %>% filter(locus == snpa | locus == snpb)
r %>% filter((SNP_A == snpa & SNP_B == snpb) | (SNP_A == snpa & SNP_B == snpb))

# Approximation of number of SNPs to retain
filtered_r %>% distinct(SNP_A) %>% nrow()
filtered_r %>% distinct(SNP_B) %>% nrow()
