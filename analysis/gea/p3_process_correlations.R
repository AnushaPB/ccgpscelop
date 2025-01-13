library(here)
library(tidyverse)
r <- read_table(here("analysis", "gea", "outputs", "snp_r2.ld"))
gea <-  read_csv(here("analysis", "gea", "outputs", "rda_sig_p01.csv"))

# Get the SNPs correlated with the RDA SNPs
rda_r <- 
  r %>% 
  filter(SNP_A %in% gea$locus | SNP_B %in% gea$locus)

# Pull out distinct SNPs
snps_r <- 
  rda_r %>%
  select(SNP_A, SNP_B) %>%
  pivot_longer(c(SNP_A, SNP_B)) %>%
  distinct(value) %>%
  pull(value)

# Combine with original SNPs
all_snps <- unique(c(gea$locus, snps_r))

# Calculate number of added SNPs
print(paste("Original number of SNPs:", nrow(gea)))
print(paste("Number of added SNPs:", length(all_snps) - nrow(gea)))

# Make into df
snp_df <- 
  tibble(locus = all_snps) %>% 
  mutate(
    # Pull out the digit in ...[digit]_[bp]_[bp] pattern
    start = as.integer(str_extract(locus, "(?<=_)[0-9]+(?=_[A-Z]+_[A-Z]+)")),
    end = start,
    scaffold = str_extract(locus, "^(Scaffold_[0-9]+__[0-9]+_contigs__length_[0-9]+|chr[0-9]+)")
  )

# Write out txt file with just ids
write.table(snp_df$locus, here("analysis", "gea", "outputs", "rda_ids.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)

# Confirm that all loci have been formatted correctly
stopifnot(all(complete.cases(snp_df)))

# Write out csv file
write_csv(snp_df, here("analysis", "gea", "outputs", "rda_ids.csv"))

# Create bed file
rda_bed <- 
  snp_df %>%
  select(scaffold, start, end) 

# Write to table with no header
write.table(rda_bed, here("analysis", "gea", "outputs", "gea.bed"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
