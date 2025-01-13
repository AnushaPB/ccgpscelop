library(here)
library(tidyverse)

rda_results <- read_csv(here("analysis", "gea", "outputs", "58-Sceloporus_RDA_outliers_full_rdadapt.csv"))
#rda_results <- read_csv(here("analysis", "gea", "outputs", "58-Sceloporus_RDA_outliers_full_Zscores.csv"))

rda_results %>% filter(grepl("Scaffold_1__1_contigs__length_4124712_470887", locus))

rda_adj <- 
  rda_results %>% 
  dplyr::rename(scaffold = scaff) %>%
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

# Write out csv file
write_csv(rda_sig, here("analysis", "gea", "outputs", "rda_sig_p01.csv"))