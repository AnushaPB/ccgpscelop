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

synonymous_by_pop <- read_table(here(outpath, "all_synonymous.frq.strat")) %>% mutate(MAF = pmin(1 - MAF, MAF))
nonsynonymous_by_pop <- read_table(here(outpath, "all_nonsynonymous.frq.strat")) %>% mutate(MAF = pmin(1 - MAF, MAF))

synonymous_by_pop %>% head(5)
nonsynonymous_by_pop %>% head(5)

breaks <- seq(0, 0.5, by = 0.01)
bin_width <- 0.01

bin_freq_by_pop <- function(df, label) {
  df %>%
    mutate(
      bin = cut(MAF,
                breaks = breaks,
                include.lowest = TRUE,
                right = FALSE,
                labels = FALSE)
    ) %>%
    group_by(CLST, bin) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(CLST) %>%
    mutate(
      prop   = count / sum(count),
      set    = label,
      x_left = breaks[bin],
      x_mid  = x_left + bin_width/2
    ) %>%
    ungroup()
}

frq_bin_counts_num %>% filter(CLST == 5, set == "Non-synonymous") %>% arrange(bin)

frq_bin_counts_num <- bind_rows(
  bin_freq_by_pop(nonsynonymous_by_pop,           "Non-synonymous"),
  bin_freq_by_pop(synonymous_by_pop,              "Synonymous")
)

source(here("general_functions.R"))
library(sf)
pop_df <- get_pops() %>% right_join(get_coords(sf = TRUE)) %>% st_as_sf()

pop_cols <- viridis::turbo(8)
names(pop_cols) <- sort(unique(pop_df$cluster))

frq_bin_by_south <- 
  frq_bin_counts_num %>% 
  mutate(
    Southern = ifelse(CLST %in% c(6, 8), "Southern", "Other")
  ) %>%
  group_by(Southern, bin, x_left, x_mid, set) %>%
    # combine pops within each region
    summarise(count = sum(count), .groups = "drop") %>% 
  group_by(Southern, set) %>% 
    # compute bin-wise proportions *within region*
    mutate(prop = count / sum(count)) %>% 
  ungroup()

pop_df <- mutate(
  pop_df,
  Southern = ifelse(cluster %in% c(6, 8), "Southern", "Other")
)
  
ca <- get_ca()


check_bins <-
  bind_rows(
    nonsynonymous_by_pop %>% mutate(set = "Non-synonymous"),
    synonymous_by_pop %>% mutate(set = "Synonymous")
  )

pdf(here(plotpath, "allele_frequency_distribution_by_pop.pdf"), width = 9, height = 4)
ggplot(check_bins) +
  geom_histogram(aes(x = MAF, fill = set), binwidth = 0.01, position = "dodge", alpha = 0.7) +
  facet_wrap(~CLST, nrow = 2, scales = "free_y")
dev.off()

plt1 <- 
  ggplot(pop_df) +
  geom_sf(data = ca) +
  geom_sf(aes(color = cluster), size = 2) +
  theme_void() +
  labs(color = "Population cluster") +
  scale_color_manual(values = pop_cols)

plt2 <-
  ggplot(frq_bin_counts_num, aes(x = x_mid, fill = set)) +
  geom_col(
    aes(y = prop),
    position = position_dodge(width = bin_width),
    width = bin_width
  ) +
  geom_text(
    data = frq_bin_counts_num %>% distinct(CLST) %>% 
      mutate(x_mid = 0.4, y = 0.6),
    aes(x = x_mid, y = y, label = CLST, color = factor(CLST)),
    inherit.aes = FALSE,
    size = 6,
    fontface = "bold"
  ) +
  scale_color_discrete(guide = "none") +  # hide text color legend
  scale_x_continuous(
    limits = c(min(breaks), max(breaks)),
    breaks = seq(0, 0.5, by = 0.05),
    labels = scales::number_format(accuracy = 0.01)
  ) +
  scale_y_continuous(expand = c(0,0)) +
  scale_color_manual(values = pop_cols) +
  guides(color = "none") +  # hide text color legend
  scale_fill_manual(values = c(
    "Non-synonymous" = "red2",
    "Synonymous" = "blue2",
    "All" = "#5f5f5f"
  )) +
  labs(
    x = "Minor Allele Frequency (MAF)",
    y = "Proportion of Variants",
    fill = "Variant set",
    title = "MAF Distribution (0.01-width bins)"
  ) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    strip.text = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) +
  facet_wrap(~ CLST, nrow = 2) +
  ylim(0, 0.7) 

pdf(here(plotpath, "allele_frequency_distribution_by_pop.pdf"), width = 15*1.5, height = 5*1.5)
cowplot::plot_grid(plt1, plt2, nrow = 1, rel_widths = c(1, 3))
dev.off()

pdf(here(plotpath, "allele_frequency_distribution_southern_vs_nonsouthern.pdf"), width = 9, height = 4)
plt1 <- 
  ggplot(pop_df) +
  geom_sf(data = ca) +
  geom_sf(aes(color = cluster), size = 2) +
  theme_void() +
  labs(color = "Population cluster") +
  scale_color_manual(values = pop_cols)

plt2 <-
  ggplot(pop_df) +
  geom_sf(data = ca) +
  geom_sf(aes(color = Southern), size = 2) +
  theme_void() +
  labs(color = "Group") 

cowplot::plot_grid(plt1, plt2, nrow = 1, align = "hv")

ggplot(frq_bin_by_south, aes(x = x_mid, y = prop, fill = set)) +
  geom_col(
    position = position_dodge(width = bin_width * 0.8),
    width = bin_width * 0.8
  ) +
  scale_x_continuous(
    limits = c(min(breaks), max(breaks)),
    breaks = seq(0, 0.5, by = 0.05),
    labels = scales::number_format(accuracy = 0.01),
    expand = c(0, 0)
  ) +
  scale_y_continuous(expand = c(0,0)) +
  labs(
    x = "Minor Allele Frequency (MAF)",
    y = "Proportion of Variants",
    fill = "Population",
    title = "MAF Distribution (0.01-width bins)"
  ) +
  theme_classic() +
  scale_color_manual(values = c("Non-synonymous" = "red2", "Synonymous" = "blue2", "All" = "#5f5f5f")) +
  facet_wrap(~ Southern) +
  theme(legend.position = "bottom")
dev.off()

