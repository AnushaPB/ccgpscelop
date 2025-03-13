library(tidyverse)
library(here)

csq_raw <- read_table("analysis/gea/outputs/all_variants.txt", col_names = FALSE)
csq <- 
  csq_raw %>%  
  rename(scaffold = X1, position = X2, csq = X3) %>%
  mutate(csq = str_extract(csq, "^[^|]*"))

write_csv(csq, here("analysis", "gea", "outputs", "csq.csv"))

nonsyn <- csq %>% filter(csq != ".", csq != "intron", csq != "synonymous")
syn <- csq %>% filter(csq == "synonymous")
exons <- csq %>% filter(csq != "intron")

write_csv(nonsyn, here("analysis", "gea", "outputs", "nonsynonymous.csv"))
write_csv(syn, here("analysis", "gea", "outputs", "synonymous.csv"))

gea_genes <- 
  read_csv(here("analysis", "gea", "outputs", "bio1ndvi_gea_gene_snp.csv")) %>%
  mutate(position = start) %>%
  left_join(exons, by = c("scaffold", "position"))

# Compare number of nonsynonymous to synonymous variants in genes
# TRUE = synonymous/other
# FALSE = nonsynonymous
table(gea_genes$csq == "synonymous")
#FALSE  TRUE 
#20,468 11,527 

# Get the non-synonymous and synonymous variants
gea_genes_syn <- gea_genes %>% filter(csq == "synonymous") %>% pull(locus)
gea_genes_nonsyn <- gea_genes %>% filter(csq != "synonymous", csq != ".") %>% pull(locus)

writeLines(gea_genes_syn, here("analysis", "gea", "outputs", "bio1ndvi_gea_gene_syn_ids.txt"))
writeLines(gea_genes_nonsyn, here("analysis", "gea", "outputs", "bio1ndvi_gea_gene_nonsyn_ids.txt"))

# Check csq of non-synonymous variants
unique(gea_genes %>% filter(csq != "synonymous", csq != ".") %>% pull(csq))
