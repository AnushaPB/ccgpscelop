library(tidyverse)
library(here)

outpath <- here("analysis", "gea", "outputs")

# Read in the csq file for all variants
csq_raw <- read_table(here(outpath, "all_variant_csq.txt"), col_names = FALSE)

# Format csq
csq <- 
  csq_raw %>%  
  rename(scaffold = X1, position = X2, csq = X3) %>%
  mutate(csq = str_extract(csq, "^[^|]*"))

# Write out the raw csq file for all variants
write_csv(csq, here(outpath, "csq.csv"))
csq <- read_csv(here(outpath, "csq.csv"))

# Identify variants with different consequences
nonsyn <- csq %>% filter(csq != ".", csq != "intron", csq != "synonymous")
syn <- csq %>% filter(csq == "synonymous")
exons <- csq %>% filter(csq != "intron")

# Write out the nonsynonymous and synonymous variants to separate files
write_csv(nonsyn, here(outpath, "nonsynonymous.csv"))
nonsyn <- read_csv(here(outpath, "nonsynonymous.csv"))

write_csv(syn, here(outpath, "synonymous.csv"))
syn <- read_csv(here(outpath, "synonymous.csv"))

# Create non-synonymous bed file
nonsyn %>%
  mutate(start = position, end = position) %>%
  select(scaffold, start, end, csq) %>%
  distinct() %>%
  write.table(here(outpath, "all_nonsynonymous.bed"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

# Read in the gene information from the GEA analysis
gea_genes <- 
  read_csv(here(outpath, "bio1ndvi_gea_gene_snp.csv")) %>%
  mutate(position = start) %>%
  left_join(exons, by = c("scaffold", "position"))

# Compare number of nonsynonymous to synonymous variants in genes
# TRUE = synonymous/other
# FALSE = nonsynonymous
table(gea_genes$csq == "synonymous")
# FALSE  TRUE 
# 31162 16792 

# Get the non-synonymous and synonymous variants
gea_genes_syn <- gea_genes %>% filter(csq == "synonymous") 
gea_genes_syn_snp <- gea_genes_syn %>% pull(locus)
gea_genes_nonsyn <- gea_genes %>% filter(csq != "synonymous", csq != ".")
gea_genes_nonsyn_snp <- gea_genes_nonsyn %>% pull(locus)

# Write out non-syn information
write_csv(gea_genes_nonsyn, here(outpath, "bio1ndvi_gea_gene_nonsyn.csv"))

# Write out the gene IDs for the synonymous and non-synonymous variants to separate files
writeLines(gea_genes_syn, here(outpath, "bio1ndvi_gea_gene_syn_ids.txt"))
writeLines(gea_genes_nonsyn, here(outpath, "bio1ndvi_gea_gene_nonsyn_ids.txt"))

# Check csq of non-synonymous variants
unique(gea_genes %>% filter(csq != "synonymous", csq != ".") %>% pull(csq))

# GET GEA SNPS NOT LINKED TO NON-SYNONYMOUS VARIANTS
r <- read_csv(here(outpath, "bio1ndvi_rda_linked_snps_info.csv"))

# Get only SNPs with linked SNPs
make_key <- function(locus){
  scaffold <- sub("_[0-9]+_[A-Z]+_[A-Z]+$", "", locus)
  position <- as.integer(sub("^.*_([0-9]+)_[A-Z]+_[A-Z]+$", "\\1", locus))
  paste0(scaffold, "_", position)
}

linked <- 
  r %>% 
  drop_na(linked_locus) %>%
  mutate(
    outlier_locus = make_key(outlier_locus),
    linked_locus = make_key(linked_locus)
  )

# Get non-synonymous CSQ 
csq <- 
  read_csv(here(outpath, "nonsynonymous.csv")) %>%
  mutate(key = paste0(scaffold, "_", position))

# Determine which SNPs are non-synonymous or linked to non-synonymous
linked_csq <-
  linked %>%
  mutate(
    nonsyn_outlier = outlier_locus %in% csq$key,
    nonsyn_linked = linked_locus %in% csq$key
  ) %>%
  mutate(
    nonsyn = nonsyn_outlier | nonsyn_linked
  )

# Pull out synonymous SNPs that are not linked to non-synonymous SNPs
linked_syn <- 
  linked_csq %>% 
  filter(!nonsyn) %>% 
  pivot_longer(c(outlier_locus, linked_locus), names_to = "type", values_to = "locus") %>%
  distinct(locus) %>%
  pull(locus)

writeLines(linked_syn, here(outpath, "bio1ndvi_gea_gene_unlinked_syn_ids.txt"))

# Get non-synonymous SNPs in genes NOT in GEA
notgeagenes <- read_csv(here(outpath, "all_genes_not_in_gea.csv"))
head(notgeagenes)

library(furrr)
# TAKES A WHILE
plan(multisession, workers = 10)
nonsyn_in_notgeagenes <- 
  notgeagenes %>%
  #select(-full_name) %>%
  filter(scaffold == "chr2") %>%
  future_pmap(\(scaffold, start, end, full_name) {
    nonsyn %>%
      filter(scaffold == scaffold, position >= start, position <= end) %>%
      mutate(full_name = full_name)
  }, .progress = TRUE) 
plan(sequential)

# Creating bed instead of text file because it is harder to recreate locus names and easier to just give intervals
notgeagenes_nonsynsnps_bed <- 
  bind_rows(nonsyn_in_notgeagenes) %>% 
  mutate(start = position, end = position) %>%
  select(scaffold, start, end, full_name) 

# Check duplicates
# notgeagenes_nonsynsnps_bed %>% 
#   filter(scaffold == "chr2") %>%
#   group_by(scaffold, start, end) %>%
#   count() %>%
#   filter(n > 1)
# notgeagenes_nonsynsnps_bed %>% filter(scaffold == "chr2", start == 1055905, end == 1055905) 
# notgeagenes %>% filter(scaffold == "chr2", 1055905 >= start, 1055905 <= end) # only one gene
# nonsyn %>% filter(scaffold == "chr2", position >= 1038546, position <= 1056579)

notgeagenes_nonsynsnps_bed_unique <- 
  notgeagenes_nonsynsnps_bed %>%
  select(scaffold, start, end) %>%
  # Remove duplicates caused by SNPs falling in multiple overlapping genes
  distinct() 

dim(distinct(notgeagenes_nonsynsnps_bed))

write.table(notgeagenes_nonsynsnps_bed_unique, here(outpath, "notgeagenes_nonsyn.bed"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
