library(tidyverse)
library(here)
library(furrr)
outpath <- here("analysis", "check_nonsyn", "outputs")

# FORMAT CSQ ------------------------------------------------------------------------------------------

# Read in the csq file for all variants
csq_raw <- read_table(here(outpath, "all_variant_csq.txt"), col_names = FALSE)
nrow(csq_raw)

# Format csq
csq <- 
  csq_raw %>%  
  rename(scaffold = X1, position = X2, csq = X3) %>%
  mutate(csq = str_extract(csq, "^[^|]*"))

csq %>% group_by(scaffold) %>% tally() %>% arrange(desc(n))

# Write out the raw csq file for all variants
write_csv(csq, here(outpath, "csq.csv"))
csq <- read_csv(here(outpath, "csq.csv"))

# Print out csq types
csq %>% group_by(csq) %>% tally() %>% arrange(desc(n))

# Identify variants with different consequences
nonsyn <- csq %>% filter(grepl("missense|stop_gained|start_lost", csq))
syn <- csq %>% filter(csq == "synonymous")
exons <- csq %>% filter(csq != "intron")

# Print number of nonsyn and syn variants
nrow(syn)
nrow(nonsyn)

# Write out the nonsynonymous and synonymous variants to separate files
write_csv(nonsyn, here(outpath, "nonsynonymous.csv"))
write_csv(syn, here(outpath, "synonymous.csv"))

# Create non-synonymous bed file
nonsyn %>%
  # Convert to 0-based start/end to join with bed files
  mutate(start = position - 1, end = position) %>%
  dplyr::select(scaffold, start, end, csq) %>%
  distinct() %>%
  mutate_at(vars(start, end), as.integer) %>%
  write_tsv(here(outpath, "all_nonsynonymous.bed"), col_names = FALSE)

# Create synonymous bed file
syn %>%
  # Convert to 0-based start/end to join with bed files
  mutate(start = position - 1, end = position) %>%
  dplyr::select(scaffold, start, end, csq) %>%
  distinct() %>%
  mutate_at(vars(start, end), as.integer) %>%
  write_tsv(here(outpath, "all_synonymous.bed"), col_names = FALSE)
