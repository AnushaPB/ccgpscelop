library(here)
library(tidyverse)
library(data.table)

# example run:
# Rscript GEA_outlier_positions.R "./outputs/58-Sceloporus_0.6.prune.out" "./outputs/58-Sceloporus_snp_correlations.ld" "./outputs/gea_genes.bed" "./outputs" "together"

#!/usr/bin/env Rscript # leave line commented
args = commandArgs(trailingOnly=TRUE)
prunedout_file = args[1]	# path to .prune.out file
corr_file = args[2] 		  # path to .ld correlations file
gea_genes_file = args[3]  # path to gea genes bed file
output_path = args[4]		  # path to save output files
output_format = args[5]		# whether to output separate bed files for pruned out SNPs vs outlier ones; options are "separate" or "together"


# Read in files -----------------------------------------------------------

prunedout <- read.table(paste0(prunedout_file), col.names = "SNP")
ldsnps <- read.table(paste0(corr_file), header = TRUE)
gea_outliers <- data.table::fread(paste0(gea_genes_file), col.names = c("CHROM", "POS"))
gea_outliers <- as.data.frame(gea_outliers)


# Process files -----------------------------------------------------------

# Unite cols in ldsnps to get CHROM_POS formatting
ldsnps <- ldsnps %>% 
  dplyr::select(-c(SNP_A, SNP_B)) %>%
  tidyr::unite("SNP_A", CHR_A:BP_A, sep = "_") %>%
  tidyr::unite("SNP_B", CHR_B:BP_B, sep = "_")

# Combine gea outliers CHROM and POS cols
gea_outliers <- gea_outliers %>%
  tidyr::unite("SNP", CHROM:POS, sep = "_")

# Split up prunedout to remove REF and ALT cols
prunedout <- prunedout %>% 
  tidyr::extract(SNP, c("CHROM", "POS", "REF", "ALT"), regex = "(.*)_([^_]+)_([^_]+)_([^_]+)$") %>% 
  dplyr::select(CHROM, POS) %>%
  tidyr::unite("SNP", CHROM:POS, sep = "_")

# Only retain correlations for outlier SNPs
corr <- ldsnps %>%
  # Only retain comparisons involving outlier SNPs
  dplyr::filter(SNP_A %in% gea_outliers$SNP | SNP_B %in% gea_outliers$SNP) %>%
  dplyr::mutate(SNP_A_category = case_when(SNP_A %in% gea_outliers$SNP ~ "outlier",
                                          SNP_A %in% prunedout$SNP ~ "prunedout"),
                SNP_B_category = case_when(SNP_B %in% gea_outliers$SNP ~ "outlier",
                                          SNP_B %in% prunedout$SNP ~ "prunedout")) %>%
  # Remove comparisons between outliers or between pruned SNPs
  dplyr::filter(SNP_A_category != SNP_B_category)

if (nrow(corr) == 0) stop("No correlations found between outliers and pruned out SNPs, stopping")

corr_a <- corr %>%
  dplyr::filter(SNP_A_category == "prunedout") %>%
  dplyr::rename(pruned_snp = SNP_A, outlier_snp = SNP_B) %>%
  dplyr::select(-c(SNP_A_category, SNP_B_category))

corr_b <- corr %>%
  dplyr::filter(SNP_B_category == "prunedout") %>%
  dplyr::rename(pruned_snp = SNP_B, outlier_snp = SNP_A) %>%
  dplyr::select(-c(SNP_A_category, SNP_B_category))

final <- bind_rows(corr_a, corr_b)


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