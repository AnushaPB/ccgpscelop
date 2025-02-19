library(here)
library(tidyverse)
library(data.table)

# example run:
# Rscript GEA_outlier_positions.R "../../data/ccgp_data/58-Sceloporus_0.6.prune.out" "../gea/outputs/snp_r2.ld" "../gea/outputs/gea_gene_ids.txt" "genes" "../gea/outputs/RDA_PCA" "together"

#!/usr/bin/env Rscript # leave line commented
args = commandArgs(trailingOnly=TRUE)
prunedout_file = args[1]	# path to .prune.out file
corr_file = args[2] 		  # path to .ld correlations file
gea_file = args[3]        # path to gea genes id file
snp_set = args[4]         # options are "rdadapt", "zscores", or "genes" (if id txt file provided)
sig = args[5]             # if `snp_set = "rdadapt"`, alpha threshold to detect outliers
outpath = args[6]         # path to output file dir
output_format = args[7]		# whether to output separate bed files for pruned out SNPs vs outlier ones; options are "separate" or "together"

### TROUBLESHOOTING - FROM APB

r <- read_table(here("analysis", "gea", "outputs", "snp_r2.ld"))
r %>% filter((SNP_A == "chr1_3757363_G_A" & SNP_B == "chr1_3757379_T_C") | (SNP_A == "chr1_3757379_T_C" & SNP_B == "chr1_3757363_G_A"))
rdasig <- read_csv(here(outpath, "RDA_bio1_ndvi", "58-Sceloporus_RDA_outliers_full_rdadapt.csv")) 
rdasig01 %>% filter(locus == "chr1_3757363_G_A" | locus == "chr1_3757379_T_C")

prunein <- read_table(here("data", "ccgp_data", "58-Sceloporus_0.6.prune.in"), col_names = "snp")
nrow(prunein)
# [1] 49306298
pruneout <- read_table(here("data", "ccgp_data", "58-Sceloporus_0.6.prune.out"), col_names = "snp")
nrow(pruneout)
# [1] 19670579
all(prunein$snp %in% pruneout$snp)
# [1] FALSE
prunein %>% filter(snp == "chr1_3757363_G_A" | snp == "chr1_3757379_T_C") # both are there
pruneout %>% filter(snp == "chr1_3757363_G_A" | snp == "chr1_3757379_T_C") # neither is there

### 

outpath = here("analysis", "gea", "outputs")
prunedout_file = here("data", "ccgp_data", "58-Sceloporus_0.6.prune.out")
corr_file = here(outpath, "snp_r2.ld")
gea_file = here(outpath, "bio1ndvi_gea_gene_ids.txt") # BIO1+NDVI
# gea_file = here(outpath, "gea_gene_ids.txt") # rasterPCs
gea_file = here(outpath, "RDA_bio1_ndvi", "58-Sceloporus_RDA_outliers_full_rdadapt.csv")
# gea_file = here(outpath, "RDA_bio1_ndvi", "58-Sceloporus_RDA_outliers_full_Zscores.csv")
# gea_file = here(outpath, "RDA_PCA", "58-Sceloporus_RDA_outliers_full_rdadapt.csv")
snp_set = "rdadapt" # "rdadapt" "zscores" "genes"
sig = 0.01

# Read in files -----------------------------------------------------------

prunedout <- read.table(paste0(prunedout_file), col.names = "SNP")
ldsnps <- read.table(paste0(corr_file), header = TRUE)
if (snp_set == "genes") {
  gea_outliers <- data.table::fread(paste0(gea_file), col.names = "SNP")
  gea_outliers <- as.data.frame(gea_outliers)
}
if (snp_set == "rdadapt") {
  gea_outliers <- data.table::fread(paste0(gea_file))
  gea_outliers <- gea_outliers %>% filter(q.values < sig) %>% rename(SNP = locus)
}

if (snp_set == "zscores") {
  gea_outliers <- data.table::fread(paste0(gea_file))
  gea_outliers <- gea_outliers %>% rename(SNP = rda_snps)
}

nrow(gea_outliers)

# Process files -----------------------------------------------------------

# Unite cols in ldsnps to get CHROM_POS formatting
# ldsnps <- ldsnps %>% 
#   dplyr::select(-c(SNP_A, SNP_B)) %>%
#   tidyr::unite("SNP_A", CHR_A:BP_A, sep = "_") %>%
#   tidyr::unite("SNP_B", CHR_B:BP_B, sep = "_")

# Combine gea outliers CHROM and POS cols
# gea_outliers <- gea_outliers %>%
#   tidyr::unite("SNP", CHROM:POS, sep = "_")

# Split up prunedout to remove REF and ALT cols
# prunedout <- prunedout %>% 
#   tidyr::extract(SNP, c("CHROM", "POS", "REF", "ALT"), regex = "(.*)_([^_]+)_([^_]+)_([^_]+)$") %>% 
#   dplyr::select(CHROM, POS) %>%
#   tidyr::unite("SNP", CHROM:POS, sep = "_")

corr <- ldsnps %>%
  # Only retain comparisons involving outlier SNPs
  dplyr::filter(SNP_A %in% gea_outliers$SNP | SNP_B %in% gea_outliers$SNP) %>%
  # Categorize SNPs based on whether they were pruned out or not
  dplyr::mutate(SNP_A_category = case_when(SNP_A %in% gea_outliers$SNP ~ "outlier",
                                          SNP_A %in% prunedout$SNP ~ "prunedout"),
                SNP_B_category = case_when(SNP_B %in% gea_outliers$SNP ~ "outlier",
                                          SNP_B %in% prunedout$SNP ~ "prunedout"))

if (nrow(corr) == 0) stop("No correlations found between outliers and pruned out SNPs, stopping")

corr_a <- corr %>%
  # Remove comparisons between outliers or between pruned SNPs
  dplyr::filter(SNP_A_category != SNP_B_category) %>%
  dplyr::filter(SNP_A_category == "prunedout") %>%
  dplyr::rename(pruned_snp = SNP_A, outlier_snp = SNP_B) %>%
  dplyr::select(-c(SNP_A_category, SNP_B_category))

corr_b <- corr %>%
  # Remove comparisons between outliers or between pruned SNPs
  dplyr::filter(SNP_A_category != SNP_B_category) %>%
  dplyr::filter(SNP_B_category == "prunedout") %>%
  dplyr::rename(pruned_snp = SNP_B, outlier_snp = SNP_A) %>%
  dplyr::select(-c(SNP_A_category, SNP_B_category))

final <- bind_rows(corr_a, corr_b)

# Were there SNPs not pruned out that exceeded r2 threshold? --------------

# Find outlier-outlier comparisons and see if r2 >= 0.6
exceed_thresh <- corr %>%
  filter(SNP_A_category == "outlier" & SNP_B_category == "outlier") %>%
  filter(R2 >= 0.6)

dist_a <- exceed_thresh %>%
  dplyr::select(SNP_A) %>%
  rename(SNP = SNP_A)
dist_b <- exceed_thresh %>%
  dplyr::select(SNP_B) %>%
  rename(SNP = SNP_B)
dat <- bind_rows(dist_a, dist_b) %>% distinct()
nrow(dat)

# Export files ------------------------------------------------------------

# Export correlation info
write_csv(final, file = paste0(output_path, "/58-Sceloporus_corr.csv"))

# Export bed files
bed_pruned <- final %>%
  dplyr::select(pruned_snp) %>%
  distinct() %>%
  tidyr::extract(pruned_snp, c("chrom", "chromStart"), regex = "(.*)_([^_]+)$") %>%
  mutate(chromEnd = chromStart)

bed_outliers <- final %>%
  dplyr::select(outlier_snp) %>%
  distinct() %>%
  tidyr::extract(outlier_snp, c("chrom", "chromStart"), regex = "(.*)_([^_]+)$") %>%
  mutate(chromEnd = chromStart)

if (output_format == "separate") {
  readr::write_tsv(bed_pruned, file = paste0(output_path, "/58-Sceloporus_gea_corr_prunedout.bed"))
  readr::write_tsv(bed_outliers, file = paste0(output_path, "/58-Sceloporus_gea_corr_outliers.bed"))
}

if (output_format == "together") {
  readr::write_tsv(bind_rows(bed_pruned, bed_outliers), file = paste0(output_path, "/58-Sceloporus_gea_corr.bed"))
}