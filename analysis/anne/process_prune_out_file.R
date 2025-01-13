library(here)
library(tidyverse)
library(data.table)

# example run:
# Rscript process_prune_out_file.R "../../data/ccgp_data/58-Sceloporus_0.6.prune.out" "./RDA_scaffolds" "./gea_genes.bed"

#!/usr/bin/env Rscript # leave line commented
args = commandArgs(trailingOnly=TRUE)
prunedout_file = args[1]	# path to .prune.out file
rda_scaffolds_file = args[2]# path to list of scaffolds run in RDA analysis
gea_genes_file = args[3]    # path to gea genes bed file

#   Original prune.out file from hgdownload `prunedout`  = 19,670,579
#       Only SNPs in relevant scaffolds `filtered`       = 19,643,155
#   Original gea_genes.bed from APB `gea_genes`          = 129,126
#       Duplicates removed                               = 121,673
#   Binding by rows pruned + gea_genes `df`              = 19,764,828
#       Duplicates removed                               = 19,763,402
#   Overlapping SNPs between `filtered` and `gea_genes`  = 1,426

# -------------------------------------------------------------------------

# Import relevant files
prunedout <- data.table::fread(paste0(prunedout_file), header = FALSE, col.names = "SNP")
rdascaffolds <- read_tsv(paste0(rda_scaffolds_file), col_names = "CHROM")
gea_genes <- data.table::fread(paste0(gea_genes_file), col.names = c("CHROM", "POS"))
gea_genes <- as.data.frame(gea_genes)
gea_genes <- gea_genes %>% distinct() # there are some duplicated SNPs, remove these

# Format prune.out file properly and extract only relevant scaffolds on which RDA was run
# Should be 20 scaffolds remaining (from 35 originally)
filtered <- prunedout %>% 
    tidyr::extract(SNP, c("CHROM", "POS", "REF", "ALT"), regex = "(.*)_([^_]+)_([^_]+)_([^_]+)$") %>% 
    dplyr::select(CHROM, POS) %>%
    dplyr::filter(CHROM %in% rdascaffolds$CHROM)
filtered$POS <- as.numeric(filtered$POS)

# Check for any overlap between gea_genes and prunedout SNPs; there should NOT be any
check_gea_genes <- gea_genes %>% unite("SNP", CHROM:POS, sep = "_")
check_filtered <- filtered %>% unite("SNP", CHROM:POS, sep = "_")
overlapping <- dplyr::inner_join(check_gea_genes, check_filtered)

df <- dplyr::bind_rows(gea_genes, filtered) %>% distinct()
write_tsv(df, "./gea_regions.txt", col_names = FALSE)
