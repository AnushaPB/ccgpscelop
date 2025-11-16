library(tidyverse)
library(here)
library(furrr)
outpath <- here("analysis", "gea", "outputs")
annotationpath <- here("analysis", "nonsyn_check", "outputs")

# FORMAT CSQ ------------------------------------------------------------------------------------------
# Read in the csq file for all variants (created using scripts in data_processing)
csq <- read_csv(here(outpath, "csq.csv"))

# Identify variants with different consequences
nonsyn <- csq %>% filter(grepl("missense|stop_gained|start_lost", csq))
syn <- csq %>% filter(csq == "synonymous")
exons <- csq %>% filter(csq != "intron")

# GET CSQ OF GEA SNPS --------------------------------------------------------------------------------

# Read in the gene information from the GEA analysis
gea_genes <- 
  read_csv(here(outpath, "bio1ndvi_gea_genes_snp.csv")) %>%
  left_join(exons, by = c("scaffold", "position"))

# Get the non-synonymous and synonymous variants
gea_genes_syn <- gea_genes %>% filter(csq == "synonymous") 
gea_genes_syn_snp <- gea_genes_syn %>% distinct(locus) %>% pull(locus)
gea_genes_nonsyn <- gea_genes %>% filter(grepl("missense|stop_gained|start_lost", csq))
gea_genes_nonsyn_snp <- gea_genes_nonsyn %>% distinct(locus) %>%  pull(locus)

# Write out non-syn information
write_csv(gea_genes_nonsyn, here(outpath, "bio1ndvi_gea_genes_nonsyn.csv"))
gea_genes_nonsyn <- read_csv(here(outpath, "bio1ndvi_gea_genes_nonsyn.csv"))
message("Number of non-synonymous SNPs in genes: ", length(gea_genes_nonsyn_snp)) #8860
message("Number of unique genes with non-synonymous SNPs: ", length(unique(gea_genes_nonsyn$full_name))) #4958

# Write out the gene IDs for the synonymous and non-synonymous variants to separate files
writeLines(gea_genes_syn_snp, here(outpath, "bio1ndvi_gea_genes_syn_ids.txt"))
writeLines(gea_genes_nonsyn_snp, here(outpath, "bio1ndvi_gea_genes_nonsyn_ids.txt"))

# Create bed file
gea_nonsyn_bed <-
  gea_genes_nonsyn %>%
  # Convert to 0-based start/end to join with bed files
  mutate(start = position - 1, end = position) %>%
  dplyr::select(scaffold, start, end) %>%
  distinct()

write.table(gea_nonsyn_bed, here(outpath, "bio1ndvi_gea_genes_nonsyn.bed"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
