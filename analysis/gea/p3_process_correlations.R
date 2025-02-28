library(here)
library(tidyverse)
r <- read_table(here("analysis", "gea", "outputs", "snp_r2.ld"))
rdasig <- read_csv(here("analysis", "gea", "outputs", "58-Sceloporus_RDA_outliers_full_rdadapt.csv")) 
rdaz <- read_csv(here("analysis", "gea", "outputs", "58-Sceloporus_RDA_outliers_full_Zscores.csv"))
#r %>% filter((SNP_A == "chr1_3757363_G_A" & SNP_B == "chr1_3757379_T_C") | (SNP_A == "chr1_3757379_T_C" & SNP_B == "chr1_3757363_G_A"))
#rdasig %>% filter(locus == "chr1_3757363_G_A" | locus == "chr1_3757379_T_C")

process_outputs <- function(input, outprefix){ 
  # Output dir
  outpath <- here("analysis", "gea", "outputs")

  # Get the SNPs correlated with the RDA SNPs
  rda_r <- 
    r %>% 
    filter(SNP_A %in% input$locus | SNP_B %in% input$locus) %>%
    select(SNP_A, SNP_B, R2)

  # Create dataframe that includes which linked SNPs are correlated with which outliers
  rda_r_out <-
    bind_rows(
      left_join(input, rda_r, by = c("locus" = "SNP_A")) %>% rename(linked_locus = SNP_B, outlier_locus = locus),
      left_join(input, rda_r, by = c("locus" = "SNP_B")) %>% rename(linked_locus = SNP_A, outlier_locus = locus)
    ) %>%
    distinct()

  write_csv(rda_r_out, here(outpath, paste0(outprefix, "_rda_linked_snps_info.csv")))

  # Pull out distinct SNPs
  snps_r <- 
    rda_r %>%
    select(SNP_A, SNP_B) %>%
    pivot_longer(c(SNP_A, SNP_B)) %>%
    distinct(value) %>%
    pull(value)

  # Combine with original SNPs
  all_snps <- unique(c(input$locus, snps_r))

  # Calculate number of added SNPs
  print(paste("Original number of SNPs:", nrow(input)))
  print(paste("Number of added SNPs:", length(all_snps) - nrow(input)))

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
  write.table(snp_df$locus, here(outpath, paste0(outprefix, "_rda_ids.txt")), quote = FALSE, row.names = FALSE, col.names = FALSE)

  # Confirm that all loci have been formatted correctly
  stopifnot(all(complete.cases(snp_df)))

  # Write out csv file with more information
  write_csv(snp_df, here(outpath,  paste0(outprefix, "_rda_ids.csv")))

  # Create bed file
  rda_bed <- 
    snp_df %>%
    select(scaffold, start, end) 

  # Write to table with no header
  write.table(rda_bed, here(outpath, paste0(outprefix, "_gea.bed")), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
}

#pca <-  read_csv(here("analysis", "gea", "outputs", "rda_sig_p01.csv"))
bio1_ndvi <-  read_csv(here("analysis", "gea", "outputs", "bio1ndvi_significant_snps.csv"))

#process_outputs(input = pca, outprefix = "pca")
process_outputs(input = bio1_ndvi, outprefix = "bio1ndvi")
 
# Check for pairs where both SNP_A and SNP_B are in rdasig$locus
# Filter pairs where both SNP_A and SNP_B are in rdasig$locus
# (e..g, pairs with R2 > 0.6 where both SNPs ended up in RDA)
filtered_r <- 
  r %>%
  filter(SNP_A %in% rdasig$locus & SNP_B %in% rdasig$locus)

# Count number of pairs
nrow(filtered_r)

# Count number of SNPs
length(unique(c(filtered_r$SNP_A, filtered_r$SNP_B)))

# Calculate summary stats on R2
mean(filtered_r$R2)
range(filtered_r$R2)


# Example:
snpa <- filtered_r %>% filter(R2 == 1) %>% slice(1) %>% pull(SNP_A)
snpb <- filtered_r %>% filter(R2 == 1) %>% slice(1) %>% pull(SNP_B)
rdasig %>% filter(locus == snpa | locus == snpb)
r %>% filter((SNP_A == snpa & SNP_B == snpb) | (SNP_A == snpa & SNP_B == snpb))
print(snpa)
print(snpb)

# Check the distances between SNPs
# Pull out the SNP position of SNP_A and SNP_B
extract_position <- function(snp) {
  as.integer(str_extract(snp, "(?<=_)[0-9]+(?=_[A-Z]+_[A-Z]+)"))
}

# Calculate the distance between SNP_A and SNP_B
snp_positions <- 
  filtered_r %>%
  mutate(
    SNP_A_pos = extract_position(SNP_A),
    SNP_B_pos = extract_position(SNP_B),
    distance = abs(SNP_A_pos - SNP_B_pos) 
  ) 

# Check the distances
mean(snp_positions$distance)
max(snp_positions$distance)
min(snp_positions$distance)
mean(snp_positions$distance < 10)
