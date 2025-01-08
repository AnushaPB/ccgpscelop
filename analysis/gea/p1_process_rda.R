library(here)
library(tidyverse)

rda_results <- read_csv(here("analysis", "gea", "outputs", "58-Sceloporus_RDA_outliers_full_rdadapt.csv"))
#rda_results <- read_csv(here("analysis", "gea", "outputs", "58-Sceloporus_RDA_outliers_full_Zscores.csv"))

rda_adj <- 
  rda_results %>% 
  rename(scaffold = scaff) %>%
  # Used a holm correction for multiple testing because it is more conservative
  mutate(p.adj = p.adjust(p.values, method = "holm")) %>%
  #filter(p.values < 0.01) %>% 
  mutate(
    # Pull out the digit in ...[digit]_[bp]_[bp] pattern
    start = as.integer(str_extract(locus, "(?<=_)[0-9]+(?=_[A-Z]+_[A-Z]+)")),
    end = start
  )

rda_sig <- rda_adj %>% filter(p.adj < 0.01)
print(paste("Number of significant loci:", nrow(rda_sig)))

# Write out txt file with just ids
write.table(rda_sig$locus, here("analysis", "gea", "outputs", "rda_ids.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)

# Confirm that all loci have been formatted correctly
stopifnot(all(complete.cases(rda_sig)))

# Write out csv file
write_csv(rda_sig, here("analysis", "gea", "outputs", "rda_sig_p01.csv"))

# Create bed file
rda_bed <- 
  rda_sig %>%
  select(scaffold, start, end) 

# Write to table with no header
write.table(rda_bed, here("analysis", "gea", "outputs", "rda_sig_p01.bed"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

# Repeat with zscores for BIO1
# NO FILTERNING NEEDED - already filtered to 3 SD
bio1_results <- read_csv(here("analysis", "gea", "outputs", "58-Sceloporus_RDA_outliers_full_Zscores.csv"))

bio1_sig <- 
  bio1_results %>% 
  rename(scaffold = scaff, locus = rda_snps) %>%
  mutate(
    # Pull out the digit in ...[digit]_[bp]_[bp] pattern
    start = as.integer(str_extract(locus, "(?<=_)[0-9]+(?=_[A-Z]+_[A-Z]+)")),
    end = start
  ) %>%
  ungroup()

# Write out txt file with just ids
write.table(bio1_sig$locus, here("analysis", "gea", "outputs", "bio1_ids.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)

# Confirm that all loci have been formatted correctly
stopifnot(all(complete.cases(bio1_sig)))

# Write out csv file
write_csv(bio1_sig, here("analysis", "gea", "outputs", "bio1_sig.csv"))

# Create bed file
bio1_bed <- 
  bio1_sig %>%
  select(scaffold, start, end) 

# Write to table with no header
write.table(bio1_bed, here("analysis", "gea", "outputs", "bio1_sig.bed"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")


