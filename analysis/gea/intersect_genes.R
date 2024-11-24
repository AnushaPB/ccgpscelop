library(tidyverse)
library(here)

# Get positions for SNPs in genes
genes <- read.table(here("analysis", "gea", "outputs", "gea_genes.bed"), sep="\t", quote="", fill=TRUE, stringsAsFactors=FALSE)
gene_pos <- 
  genes %>% 
  dplyr::select(V1, V2, V3, V12) %>% 
  dplyr::rename(scaffold = V1, start = V2, end = V3, full_name = V12) %>%
  # Remove duplicates (if you don't you will get an error later with the joins)
  distinct()

# Get RDA positions which have SNP names
rda <- 
  read_csv(here("analysis", "gea", "outputs", "rda_sig_p01.csv")) %>% 
  dplyr::select(scaffold, start, end, locus, p.adj)

# Intersect the two to get SNP names for the SNPs in genes
gene_snp <- left_join(gene_pos, rda, by = c("scaffold", "start", "end"))

# Write out just the locus names
write.table(gene_snp$locus, here("analysis", "gea", "outputs", "gene_ids.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)

write.csv(gene_snp, here("analysis", "gea", "outputs", "gene_snp.csv"), row.names = FALSE)
