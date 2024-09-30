library(here)
library(tidyverse)

#!/usr/bin/env Rscript # leave line commented
args = commandArgs(trailingOnly=TRUE)
prunedout_file = args[1]	# path to .prune.out file
corr_file = args[2]		    # path to .ld correlations file
gea_path = args[3]		    # path to concatenated GEA results
sig = args[4]             # alpha threshold for rdadapt p-values
output_format = args[5]		# whether to output separate bed files for pruned out SNPs vs outlier ones; options are "separate" or "together"
output_path = args[6]		  # output path


# Read in files -----------------------------------------------------------

prunedout <- read.table(paste0(prunedout_file), col.names = "SNP")
ldsnps <- read.table(paste0(corr_file), header = TRUE)
zscores <- read_csv(paste0(gea_path, "/58-Sceloporus_Zscores.csv"),
                    col_names = c("rda_snps", "axis", "loading"))
rdadapt <- read_csv(paste0(gea_path, "/58-Sceloporus_rdadapt.csv"),
                    col_names = c("p.values", "q.values", "rda_snps")) %>%
  filter(p.values <= sig)


# Process files -----------------------------------------------------------

# Only retain correlations for outlier SNPs
process_corr <- function(dat, rda_snps, prunedout) {
  corr <- dat %>%
    dplyr::filter(SNP_A %in% rda_snps |
                    SNP_B %in% rda_snps) %>%
    dplyr::mutate(SNP_A_category = case_when(SNP_A %in% rda_snps ~ "outlier",
                                             SNP_A %in% prunedout$SNP ~ "prunedout"),
                  SNP_B_category = case_when(SNP_B %in% rda_snps ~ "outlier",
                                             SNP_B %in% prunedout$SNP ~ "prunedout")) %>%
    dplyr::filter(SNP_A_category != SNP_B_category)

  if (nrow(corr) == 0) message("No correlations found between outliers and pruned out SNPs, stopping")

  corr_a <- corr %>%
    dplyr::filter(SNP_A_category == "prunedout") %>%
    dplyr::rename(pruned_chrom = CHR_A,
                  pruned_pos = BP_A,
                  pruned_snp = SNP_A,
                  outlier_chrom = CHR_B,
                  outlier_pos = BP_B,
                  outlier_snp = SNP_B) %>%
    dplyr::select(-c(SNP_A_category, SNP_B_category))

  corr_b <- corr %>%
    dplyr::filter(SNP_B_category == "prunedout") %>%
    dplyr::rename(pruned_chrom = CHR_B,
                  pruned_pos = BP_B,
                  pruned_snp = SNP_B,
                  outlier_chrom = CHR_A,
                  outlier_pos = BP_A,
                  outlier_snp = SNP_A) %>%
    dplyr::select(-c(SNP_A_category, SNP_B_category))

  final <- bind_rows(corr_a, corr_b)

  return(final)
}

# Z-scores
corr_z <- process_corr(dat = ldsnps, rda_snps = zscores$rda_snps, prunedout = prunedout)

# RDadapt SNPs
corr_p <- process_corr(dat = ldsnps, rda_snps = rdadapt$rda_snps, prunedout = prunedout)

if (nrow(corr_z) == 0) corr_z <- NULL
if (nrow(corr_p) == 0) corr_p <- NULL


# Export files ------------------------------------------------------------

# Export correlation info
if (!is.null(corr_z)) {
  write_csv(corr_z, file = paste0(output_path, "/", species, "_corr_Zscores.csv"))
}
if (!is.null(corr_p)) {
  write_csv(corr_p, file = paste0(output_path, "/", species, "_corr_rdadapt.csv"))
}

# Export bed files
bed_helper <- function(corr) {
  bed_pruned <- corr %>%
    dplyr::select(pruned_chrom, pruned_pos, pruned_snp) %>%
    distinct() %>%
    rename(chrom = pruned_chrom,
           chromStart = pruned_pos) %>%
    mutate(chromEnd = chromStart) %>%
    dplyr::select(-pruned_snp)
  return(bed_pruned)
}

if (!is.null(corr_z)) {
  bed_z <- bed_helper(corr_z)
}
if (!is.null(corr_p)) {
  bed_p <- bed_helper(corr_p)
}

bed_rdadapt <- rdadapt %>%
  tidyr::separate(rda_snps, sep = "_", into = c("chrom", "chromStart", "ref", "alt")) %>%
  dplyr::mutate(chromEnd = as.numeric(chromStart),
                chromStart = as.numeric(chromStart)) %>%
  dplyr::select(-c(ref, alt, q.values, p.values))

bed_zscores <- zscores %>%
  tidyr::separate(rda_snps, sep = "_", into = c("chrom", "chromStart", "ref", "alt")) %>%
  dplyr::mutate(chromEnd = as.numeric(chromStart),
                chromStart = as.numeric(chromStart)) %>%
  dplyr::select(-c(ref, alt, axis, loading))

if (output_format == "separate") {
  if (!is.null(bed_p)) write_tsv(bed_p, file = paste0(output_path, "/", species, "_gea_prunedout_rdadapt.bed"))
  if (!is.null(bed_z)) write_tsv(bed_z, file = paste0(output_path, "/", species, "_gea_prunedout_Zscores.bed"))
  write_tsv(bed_rdadapt, file = paste0(output_path, "/", species, "_gea_rda_rdadapt.bed"))
  write_tsv(bed_zscores, file = paste0(output_path, "/", species, "_gea_rda_Zscores.bed"))
}

# TODO didn't account for if bed_p or bed_z are NULL below
if (output_format == "together") {
  write_tsv(bind_rows(bed_p, bed_rdadapt), file = paste0(output_path, "/", species, "_gea_rda_rdadapt.bed"))
  write_tsv(bind_rows(bed_z, bed_zscores), file = paste0(output_path, "/", species, "_gea_rda_Zscores.bed"))
}
