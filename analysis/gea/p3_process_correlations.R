library(here)
library(tidyverse)

process_outputs <- function(input, outprefix, r) {
  options(scipen = 999)  # disable scientific notation for writes
  
  # Output dir
  outpath <- here("analysis", "gea", "outputs")
  dir.create(outpath, showWarnings = FALSE, recursive = TRUE)

  # --- 1) Get pairs in LD with outlier SNPs (rda-linked) ---
  rda_r <- r %>%
    filter(SNP_A %in% input$locus | SNP_B %in% input$locus) %>%
    select(SNP_A, SNP_B, R2)

  rda_r_out <- bind_rows(
    left_join(input, rda_r, by = c("locus" = "SNP_A")) %>%
      rename(linked_locus = SNP_B, outlier_locus = locus),
    left_join(input, rda_r, by = c("locus" = "SNP_B")) %>%
      rename(linked_locus = SNP_A, outlier_locus = locus)
  ) %>%
    distinct()

  write_csv(rda_r_out, file.path(outpath, paste0(outprefix, "_rda_linked_snps_info.csv")))

  # --- 2) Collect the union of outlier SNPs + linked SNPs ---
  snps_r <- rda_r %>%
    select(SNP_A, SNP_B) %>%
    pivot_longer(c(SNP_A, SNP_B), values_to = "locus") %>%
    distinct(locus) %>%
    pull(locus)

  all_snps <- unique(c(input$locus, snps_r))
  all_snps <- all_snps[!is.na(all_snps)]

  message("Original number of SNPs (unique): ", length(unique(input$locus)))
  message("Number of added SNPs (unique): ", length(all_snps) - length(unique(input$locus)))

  # --- 3) Parse locus into scaffold + position robustly ---
  # Assumes locus format contains "..._<position>_<allele>_<allele>"
  # Works with scaffold names that themselves contain underscores (e.g., SCAF_13, Scaffold_99__...).
  snp_df <- tibble(locus = all_snps) %>%
    mutate(
      position = str_extract(locus, "(?<=_)\\d+(?=_[A-Za-z]+_[A-Za-z]+)"),
      position = suppressWarnings(as.integer(position)),
      scaffold = str_remove(locus, "_\\d+_[A-Za-z]+_[A-Za-z]+$")  # keep everything before _<pos>_<allele>_<allele>
    )

  # sanity check with helpful error
  bad_rows <- which(!complete.cases(snp_df))
  if (length(bad_rows) > 0) {
    stop(
      "Failed to parse locus → (scaffold, position) for ",
      length(bad_rows), " record(s). Example:\n",
      paste0(utils::capture.output(print(head(snp_df[bad_rows, ], 10))), collapse = "\n"),
      "\nCheck locus naming; expected pattern: <scaffold>_<pos>_<allele>_<allele>"
    )
  }

  # --- 4) Write ID lists ---
  write.table(
    snp_df$locus,
    file = file.path(outpath, paste0(outprefix, "_rda_ids.txt")),
    quote = FALSE, row.names = FALSE, col.names = FALSE
  )
  write_csv(snp_df, file.path(outpath, paste0(outprefix, "_rda_ids.csv")))

  # --- 5) Make BED (0-based, half-open, width 1) ---
  rda_bed <- snp_df %>%
    transmute(
      scaffold = scaffold,
      start = as.integer(position - 1L),
      end   = as.integer(position)
    )

  # final safety: ensure no NA and start < end
  stopifnot(all(complete.cases(rda_bed)))
  stopifnot(all(rda_bed$start < rda_bed$end))

  write.table(
    rda_bed,
    file = file.path(outpath, paste0(outprefix, "_gea.bed")),
    quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t"
  )
}

# Get SNP correlations
r <- read_table(here("analysis", "gea", "outputs", "snp_r2.ld"))

#pca <-  read_csv(here("analysis", "gea", "outputs", "rda_sig_p01.csv"))
bio1_ndvi <-  read_csv(here("analysis", "gea", "outputs", "bio1ndvi_significant_snps_unlinked.csv"))

#process_outputs(input = pca, outprefix = "pca")
process_outputs(input = bio1_ndvi, outprefix = "bio1ndvi", r = r)
# [1] "Original number of SNPs (unique): 1449874"                               
# [1] "Number of added SNPs (unique): 421577"
 
# Check for pairs where both SNP_A and SNP_B are in rdasig$locus
# Filter pairs where both SNP_A and SNP_B are in rdasig$locus
# (e..g, pairs with R2 > 0.6 where both SNPs ended up in RDA)
filtered_r <- 
  r %>%
  filter(SNP_A %in% bio1_ndvi$locus & SNP_B %in% bio1_ndvi$locus)

# Count number of pairs
nrow(filtered_r) # 4869 # COME BACK TO THIS

# Count number of SNPs
length(unique(c(filtered_r$SNP_A, filtered_r$SNP_B))) #9647

# Calculate summary stats on R2
mean(filtered_r$R2)
range(filtered_r$R2)

# Example:
snpa <- filtered_r %>% filter(R2 > 0.6) %>% slice(1) %>% pull(SNP_A)
snpb <- filtered_r %>% filter(R2 > 0.6) %>% slice(1) %>% pull(SNP_B)
bio1_ndvi %>% filter(locus == snpa | locus == snpb)
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

