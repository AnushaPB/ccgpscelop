library(tidyverse)
library(here)
outpath <- here("analysis", "gea", "outputs")
genes <- read_csv(here(outpath, "bio1ndvi_gea_gene_snp.csv"))

cortest <- 
  read_csv(here(outpath, "RDA_bio1_ndvi", "58-Sceloporus_RDA_cortest_full.csv")) %>% 
  dplyr::rename(locus = snp) %>%
  filter(locus %in% genes$locus) %>%
  filter(outlier_method == "p") %>%
  dplyr::select(r, p, var, locus)

# Note: some loci in genes will not be in cortest because they are the linked loci
# Conversely, some loci in cortest will not be in genes because they are not in genes
rda_r <- 
  cortest %>%
  # Need to filter before grouping
  filter(p < 0.05) %>%
  dplyr::select(r, var, locus) %>%
  pivot_wider(names_from = var, values_from = r) %>%
  mutate(group = case_when(
    !is.na(`CHELSA_bio1_1981.2010_V.2.1`) & is.na(NDVI) ~ "bio1",
    is.na(`CHELSA_bio1_1981.2010_V.2.1`) & !is.na(NDVI) ~ "ndvi",
    !is.na(`CHELSA_bio1_1981.2010_V.2.1`) & !is.na(NDVI) ~ "both"
  ))

rda_r %>% group_by(group) %>% count()

# Bring in the info about linked loci
linked_info <- 
  read_csv(here(outpath, "bio1ndvi_rda_linked_snps_info.csv")) %>%
  dplyr::select(outlier_locus, linked_locus) %>%
  # Note: some outlier loci do not have linked loci
  drop_na(linked_locus)

joined <- 
  left_join(linked_info, rda_r, by = c("outlier_locus" = "locus")) %>%
  # Filter to loci in genes
  dplyr::filter(linked_locus %in% genes$locus & outlier_locus %in% genes$locus) %>%
  # Drop outlier_locus column and rename linked_locus to locus and bind to rda_r
  dplyr::select(-outlier_locus) %>%
  dplyr::rename(locus = linked_locus) %>%
  bind_rows(rda_r) 

write_csv(joined, here(outpath, "bio1ndvi_genes_cortest.csv"))

# FIGURING OUT MYSTERY

outpath <- here("analysis", "gea", "outputs")
rdasig <- read_csv(here(outpath, "RDA_bio1_ndvi", "58-Sceloporus_RDA_outliers_full_rdadapt.csv")) 
rdasig05 <- rdasig %>% filter(p.values < 0.05)
cortest <- read_csv(here(outpath, "RDA_bio1_ndvi", "58-Sceloporus_RDA_cortest_full.csv"))
all(rdasig$locus %in% cortest$snp) 
# FALSE
all(rdasig05$locus %in% cortest$snp) 
# FALSE

nrow(rdasig)
#[1] 50071003
nrow(cortest)
#[1] 5405596

diff_loci <- setdiff(rdasig$locus, cortest$snp) 
length(diff_loci) # 48,558,133
diff_loci[1]
