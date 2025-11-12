library(tidyverse)
library(here)
outpath <- here("analysis", "check_nonsyn", "outputs")
plotpath <- here("analysis", "check_nonsyn", "plots")

# pmin(1 - MAF, MAF) to get MAF from allele frequency (i.e., make sure your MAF is for the minor allele)
allele_frq <- read_table(here(outpath, "58-Sceloporus.frq")) %>% mutate(MAF = pmin(1 - MAF, MAF))
nonsyn_frq <- read_table(here(outpath, "all_nonsynonymous.frq")) %>% filter(CHR %in% allele_frq$CHR) %>% mutate(MAF = pmin(1 - MAF, MAF))
syn_frq <- read_table(here(outpath, "all_synonymous.frq")) %>% filter(CHR %in% allele_frq$CHR) %>% mutate(MAF = pmin(1 - MAF, MAF))
# You can ignore warning for SCAF since we only use chromosomes

# Calculate allele frequency distribution in 0.01-width bins from 0 to 0.5
breaks <- seq(0, 0.5, by = 0.01)
bin_width <- 0.01

bin_freq <- function(df, label) {
  idx <- cut(df$MAF, breaks = breaks, include.lowest = TRUE, right = FALSE, labels = FALSE)
  tibble(bin = idx) %>%
    count(bin, name = "count") %>%
    mutate(
      prop    = count / sum(count),
      set     = label,
      x_left  = breaks[bin],
      x_mid   = x_left + bin_width/2
    )
}

frq_bin_counts_num <- bind_rows(
  bin_freq(allele_frq_subsample, "All"),
  bin_freq(nonsyn_frq,           "Non-synonymous"),
  bin_freq(syn_frq,              "Synonymous")
)

frq_bin_counts_num %>% 
  filter(set == "Non-synonymous") %>%
  arrange(count)

frq_bin_counts_num %>% 
  filter(set == "All") %>%
  arrange(count)

pdf(here(plotpath, "allele_frequency_distribution.pdf"), width=6, height=4)

ggplot(frq_bin_counts_num, aes(x = x_mid, y = prop, fill = set)) +
  geom_col(
    position = position_dodge(width = bin_width * 0.8),
    width = bin_width * 0.8,
    alpha = 0.7
  ) +
  scale_x_continuous(
    limits = c(min(breaks), max(breaks)),
    breaks = seq(0, 0.5, by = 0.05),
    labels = scales::number_format(accuracy = 0.01),
    expand = c(0, 0)
  ) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = c("Non-synonymous" = "red2", "Synonymous" = "blue2", "All" = "#5f5f5f")) +
  labs(
    x = "Minor Allele Frequency (MAF)",
    y = "Proportion of Variants",
    fill = "Variant set",
    title = "MAF Distribution (0.01-width bins)"
  ) +
  theme_classic() +
  theme(legend.position = "bottom")
dev.off()

# Define 0.05-width bins from 0 to 0.5
breaks <- seq(0, 0.5, by = 0.01)
labels <- paste0(sprintf("%.2f", head(breaks, -1)), "–", sprintf("%.2f", tail(breaks, -1)))

bin_freq <- function(df, label) {
  df %>%
    mutate(bin = cut(MAF,
                     breaks = breaks,
                     include.lowest = TRUE,
                     right = FALSE,
                     labels = labels)) %>%
    count(bin, name = "count") %>%
    mutate(prop = count / sum(count),
           set  = label)
}

binned_all    <- bin_freq(allele_frq, "all_sites")
binned_nonsyn <- bin_freq(nonsyn_frq,  "nonsyn")
binned_syn    <- bin_freq(syn_frq,  "syn")

binned_all %>% head(5)
binned_nonsyn %>% head(5)
binned_syn %>% head(5)

# Extract variants in each 0.01 bin from 0 to 0.5
variants_by_bin <- 
  allele_frq %>%                
  mutate(bin = cut(MAF,
                   breaks = breaks,
                   include.lowest = TRUE,
                   right = FALSE,
                   labels = labels)) %>%
  group_by(bin) %>%
  nest()                                   # makes a list-column with the variants in each bin

tail(variants_by_bin)

walk2(
  .x = variants_by_bin$data,
  .y = variants_by_bin$bin,
  ~ write_tsv(
      .x %>% select(SNP),   # PLINK-friendly: SNP IDs only
      here(outpath, paste0("variants_MAF_", .y, ".tsv")),
      col_names = FALSE     # PLINK extract files should NOT have a header
    )
)

nonsyn_variants_by_bin <- 
  nonsyn_frq %>%
  mutate(bin = cut(MAF,
                   breaks = breaks,
                   include.lowest = TRUE,
                   right = FALSE,
                   labels = labels)) %>%
  group_by(bin) %>%
  nest()                                   # makes a list-column with the variants in each bin

walk2(
  .x = nonsyn_variants_by_bin$data,
  .y = nonsyn_variants_by_bin$bin,
  ~ write_tsv(
      .x %>% select(SNP),   # PLINK-friendly: SNP IDs only
      here(outpath, paste0("nonsyn_variants_MAF_", .y, ".tsv")),
      col_names = FALSE     # PLINK extract files should NOT have a header
    )
)

syn_variants_by_bin <- 
  syn_frq %>%
  mutate(bin = cut(MAF,
                   breaks = breaks,
                   include.lowest = TRUE,
                   right = FALSE,
                   labels = labels)) %>%
  group_by(bin) %>%
  nest()                                   # makes a list-column with the variants in each bin  

walk2(
  .x = syn_variants_by_bin$data,
  .y = syn_variants_by_bin$bin,
  ~ write_tsv(
      .x %>% select(SNP),   # PLINK-friendly: SNP IDs only
      here(outpath, paste0("syn_variants_MAF_", .y, ".tsv")),
      col_names = FALSE     # PLINK extract files should NOT have a header
    )
)

head(nonsyn_variants_by_bin)
head(syn_variants_by_bin)
tail(nonsyn_variants_by_bin)
tail(syn_variants_by_bin)
