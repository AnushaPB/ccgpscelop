# Export complete list of genes with nonsynonymous GEA variants for Dryad
library(tidyverse)
library(here)
source(here("analysis", "gea", "functions_selection_stats.R"))

# Get genes containing nonsynonymous GEA variants
genes_nonsyn <- get_gea_genes(nonsyn = TRUE)

# Create complete gene list
gene_list <- genes_nonsyn %>%
  filter(!grepl("egapxtmp", gene_name)) %>%
  distinct(gene_name, full_name) %>%
  arrange(gene_name)

# Export for Dryad
write_csv(gene_list, here("analysis", "gea", "outputs", "nonsynonymous_gea_genes.csv"))

cat("Exported", nrow(gene_list), "distinct genes\n")
