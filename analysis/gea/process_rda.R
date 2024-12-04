library(here)
library(tidyverse)

folders <- list.files(here("outputs", "RDA"))
names(folders) <- folders
files <- map_chr(folders, ~list.files(here("outputs", "RDA", .x), full.names = TRUE, pattern = "rdadapt.csv"))
rda_results <- map(files, ~read.csv(.x), .progress = TRUE) %>% bind_rows(.id = "scaffold")
rda_sig <- 
  rda_results %>% 
  # Used a holm correction for multiple testing because it is more conservative
  mutate(p.adj = p.adjust(p.values, method = "holm")) %>%
  filter(p.values < 0.01)

# Write out txt file with just ids
write.table(rda_sig$locus, here("analysis", "gea", "outputs", "rda_ids.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)

# Extract numeric coordinates from locus
rda_sig <- 
  rda_sig %>% 
  mutate(
    start = as.integer(str_extract(locus, "(?<=_)[0-9]+")),
    end = as.integer(str_extract(locus, "(?<=_)[0-9]+"))
  )

# Write out csv file
write_csv(rda_sig, here("analysis", "gea", "outputs", "rda_sig_p01.csv"))

# Create bed file
rda_bed <- 
  rda_sig %>%
  select(scaffold, start, end) 

# Write to table with no header
write.table(rda_bed, here("analysis", "gea", "outputs", "rda_sig_p01.bed"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

library(vcfR)
library(here)
vcf <- read.vcfR(here("data", "ccgp_data", "58-Sceloporus_annotated_pruned_0.6.vcf.gz"))
