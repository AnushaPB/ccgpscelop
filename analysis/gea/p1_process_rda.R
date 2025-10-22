library(here)
library(tidyverse)

# BIO1 + NDVI
chrs <- paste0("chr", 1:10)
names(chrs) <- chrs
rda_results <- 
  map(chrs, ~read_csv(here("analysis", "gea", "outputs", "RDA_results", .x, "58-Sceloporus_RDA_outliers_full_Zscores.csv"))) %>%
  bind_rows(.id = "chr")

rda_results <- 
  rda_results %>% 
  dplyr::rename(scaffold = chr, locus = rda_snps) %>%
  mutate(
    # Pull out the digit in ...[digit]_[bp]_[bp] pattern
    position = as.integer(str_extract(locus, "(?<=_)[0-9]+(?=_[A-Z]+_[A-Z]+)")),
  ) 

# IMPORTANT: remove duplicate SNPs (e.g., SNPs associated with more than one axis)
rda_snps <-
  rda_results %>%
  distinct(scaffold, locus, position)

print(paste("Number of significant loci:", nrow(rda_snps)))
# "Number of significant loci: 1602968"

# Write out csv file of SNPs
write_csv(rda_snps, here("analysis", "gea", "outputs", "bio1ndvi_significant_snps_unlinked.csv"))

# Write out csv file of Z-scores
write_csv(rda_results, here("analysis", "gea", "outputs", "bio1ndvi_zscores.csv"))
