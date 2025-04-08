library(here)
library(tidyverse)

# BIO1 + NDVI
rda_results <- read_csv(here("analysis", "gea", "outputs", "58-Sceloporus_RDA_outliers_full_Zscores.csv"))

rda_results <- 
  rda_results %>% 
  dplyr::select(-path) %>%
  dplyr::rename(scaffold = scaff, locus = rda_snps) %>%
  mutate(
    # Pull out the digit in ...[digit]_[bp]_[bp] pattern
    start = as.integer(str_extract(locus, "(?<=_)[0-9]+(?=_[A-Z]+_[A-Z]+)")),
    end = start
  ) %>%
  # IMPORTANT: remove duplicate SNPs (e.g., SNPs associated with more than one axis)
  distinct(locus, axis, loading)

print(paste("Number of significant loci:", nrow(rda_results)))
#  "Number of significant loci: 1602968"

# Write out csv file
write_csv(rda_results, here("analysis", "gea", "outputs", "bio1ndvi_significant_snps.csv"))

# # Load rdadapt results
# rda_results <- read_csv(here("analysis", "gea", "outputs", "58-Sceloporus_RDA_outliers_full_rdadapt.csv"))

# rda_adj <- 
#   rda_results %>% 
#   dplyr::rename(scaffold = scaff) %>%
#   # Used a holm correction for multiple testing because it is more conservative
#   mutate(p.adj = p.adjust(p.values, method = "holm")) %>%
#   mutate(
#     # Pull out the digit in ...[digit]_[bp]_[bp] pattern
#     start = as.integer(str_extract(locus, "(?<=_)[0-9]+(?=_[A-Z]+_[A-Z]+)")),
#     end = start
#   )

# rda_sig <- rda_adj %>% filter(p.adj < 0.01)
# print(paste("Number of significant loci:", nrow(rda_sig)))

# # Write out csv file
# write_csv(rda_sig, here("analysis", "gea", "outputs", "rda_sig_p01.csv"))