library(tidyverse)
library(here)
source(here("general_functions.R"))
source(here("analysis", "gea", "functions_selection_stats.R"))
outpath <- here("analysis", "gea", "outputs")

# Get all non-synonymous variants
all_nonsyn <- read_csv(here(outpath, "nonsynonymous.csv"))

# Get non-synonymous GEA variants
genes_nonsyn <- read_csv(here(outpath, "bio1ndvi_gea_genes_nonsyn.csv"))

# Get all annotated genes
all_genes <-read_csv(here("analysis", "gea", "outputs", "all_genes_list.csv"))

# Filter to only include genes with non-synonymous GEA variants
genes_org <- 
  all_genes %>% 
  # TODO: this shouldn't really need to be two steps
  filter(full_name %in% genes_nonsyn$full_name) 

print(paste("Number of  Non-synonymous GEA SNPs in genes:", nrow(genes_nonsyn))) # 8,860
print(paste("Number of unique genes with Non-synonymous GEA SNPs:", length(unique(genes_nonsyn$full_name)))) # 4,958
print(paste("Number of unique Non-synonymous GEA gene names:", nrow(genes_org %>% distinct(gene_name)))) # 3,522
print(paste("Number of all genes:", nrow(all_genes))) # 18,670
print(paste("Number of all unique gene names:", nrow(all_genes %>% distinct(gene_name)))) # 13,136
print(paste("Number of non-synonymous SNPs:", nrow(all_nonsyn))) # 216,156

# PLOTTING RDA
# IMPORTANT: NEED TO DO BY SCAFFOLD!!!!!!
chrs <- paste0("chr", 1:10)
names(chrs) <- chrs
rda_z <- 
  map(chrs, ~read_csv(here("analysis", "gea", "outputs", "RDA_results", .x, "58-Sceloporus_unscaledloadings.csv"))) %>% bind_rows(.id = "chr")

# Filter to just locus scores
rda_z <- 
  rda_z %>% 
  filter(score == "species") %>% 
  rename(locus = label) %>%
  # RESCALING FOR PLOTTING
  group_by(chr) %>%
  mutate(across(starts_with("RDA"),  ~ as.vector(scale(.))))

rda_varloadings <- 
  map(chrs, ~ read_csv(here("analysis", "gea", "outputs", "RDA_results", .x, "58-Sceloporus_biplot.csv"))) %>%
  bind_rows(.id = "chr") %>%
  mutate(var = make_pretty_names(var)) %>%
  mutate(var = gsub("Contemporary t", "T", var)) %>%
  mutate(chr = factor(chr, levels = paste0("chr", 1:10))) 

# Get RDA loadings for significant loci
# TO DO FIGURE OUT DIFFERENCE BETWEEN THESE TWO:
gea_sig1 <- read_csv(here(outpath, "bio1ndvi_rda_ids.csv"))
gea_sig2 <- read_csv(here(outpath, "bio1ndvi_significant_snps_unlinked.csv"))

gea <- gea_sig2
rda_sig_z <- 
  rda_z %>%
  filter(locus %in% gea$locus) 

# Get loadings for SNPs of interest
# This loads SNPs that are non-synonymous in genes of interest, but gives only original outlier loci + any genes associated with that outlier loci AND any linked snps
genes_nonsyn <- get_cor_snps_info()
rda_genes <-
  genes_nonsyn %>%
  right_join(rda_sig_z) 

# Regex pattern for fast binary match
goi <- get_goi_names()
goi <- goi[goi != "Heat shock 70"] # Remove general HSP70
pattern <- str_c(goi, collapse = "|")
goi_df <- data.frame(goi_name = names(goi), goi = goi, stringsAsFactors = FALSE)
gg_df <- 
  rda_genes %>%
  mutate(
    goi = ifelse(str_detect(full_name, pattern), str_extract(full_name, pattern), NA)
  ) %>%
  left_join(goi_df) %>%
  drop_na(RDA1, RDA2)

goi_gg_df <- gg_df %>% filter(!is.na(goi_name))

# Calculate 3 SD cutoffs
cutoffs <-
  rda_z %>%
  group_by(chr) %>%
  summarise(
    mean_RDA1 = mean(RDA1, na.rm = TRUE),
    mean_RDA2 = mean(RDA2, na.rm = TRUE),
    sd_RDA1 = sd(RDA1, na.rm = TRUE),
    sd_RDA2 = sd(RDA2, na.rm = TRUE)
  ) %>%
  mutate(
    RDA1_3sd_max = mean_RDA1 + (3 * sd_RDA1),
    RDA2_3sd_max = mean_RDA2 + (3 * sd_RDA2),
    RDA1_3sd_min = mean_RDA1 - (3 * sd_RDA1),
    RDA2_3sd_min = mean_RDA2 - (3 * sd_RDA2)
  )

# Top 5 RDA loadings for each chromosome and each RDA axis
top_rda1 <- 
  rda_z %>% 
  group_by(chr) %>% 
  slice_max(RDA1, n = 100, with_ties = FALSE) %>% 
  ungroup()
top_rda2 <- 
  rda_z %>% 
  group_by(chr) %>% 
  slice_max(RDA2, n = 100, with_ties = FALSE) %>% 
  ungroup()
top_rda <- bind_rows(top_rda1, top_rda2) %>% distinct()
genes_nonsyn %>% filter(locus %in% top_rda$locus) %>% drop_na(gene_name) %>% dplyr::select(gene_name) %>% print(n = nrow(.)) # 1


# Top 5 for plotting
top_rda1 <- 
  rda_z %>% 
  group_by(chr) %>% 
  slice_max(RDA1, n = 5, with_ties = FALSE) %>% 
  ungroup()
top_rda2 <- 
  rda_z %>% 
  group_by(chr) %>% 
  slice_max(RDA2, n = 5, with_ties = FALSE) %>% 
  ungroup()
top_rda <- bind_rows(top_rda1, top_rda2) %>% distinct()

rda_varloadings_mean <- 
  rda_varloadings %>%
  group_by(var) %>%
  summarise(
    RDA1 = mean(RDA1, na.rm = TRUE),
    RDA2 = mean(RDA2, na.rm = TRUE)
  )

plotpath <- here("analysis", "gea", "plots")
pdf(here(plotpath, "goi_rda.pdf"), width = 6, height = 5)
zscale = 8

ggplot() +
  #geom_point(data = rda_sample, aes(x = RDA1, y = RDA2), col = "gray", alpha = 0.5) +
  geom_hex(data = rda_z, aes(x = RDA1, y = RDA2), binwidth = 0.2, alpha = 1) +
  geom_text(data = rda_varloadings_mean, aes(x = RDA1 * zscale, y = RDA2 * zscale, label = var), size = 4, vjust = 1, hjust = 1) +
  geom_segment(data = rda_varloadings_mean, aes(x = 0, y = 0, xend = RDA1 * zscale, yend = RDA2 * zscale),
               arrow = arrow(length = unit(0.2, "cm"))) +
  scale_color_identity() +
  scale_fill_gradientn(
    trans = "log10",
    colours = gray.colors(100, start = 1, end = 0)
  ) +
  geom_point(data = goi_gg_df, aes(x = RDA1, y = RDA2), col = "coral") +
  ggrepel::geom_text_repel(
    data = goi_gg_df,
    aes(x = RDA1, y = RDA2, label = goi_name),
    size = 4, col = "coral", max.overlaps = Inf
  ) +
  labs(col = "Genes of interest", fill = "SNP\ncount") +
  theme_void() 

ggplot() +
  #geom_point(data = rda_sample, aes(x = RDA1, y = RDA2), col = "gray", alpha = 0.5) +
  geom_hex(data = rda_z, aes(x = RDA1, y = RDA2), binwidth = 0.2, alpha = 1) +
  geom_text(data = rda_varloadings_mean, aes(x = RDA1 * zscale, y = RDA2 * zscale, label = var), size = 4, vjust = 1, hjust = 1) +
  geom_segment(data = rda_varloadings_mean, aes(x = 0, y = 0, xend = RDA1 * zscale, yend = RDA2 * zscale),
               arrow = arrow(length = unit(0.2, "cm"))) +
  scale_color_identity() +
  scale_fill_gradientn(
    colours = gray.colors(100, start = 1, end = 0)
  ) +
  geom_point(data = goi_gg_df, aes(x = RDA1, y = RDA2), col = "coral") +
  ggrepel::geom_text_repel(
    data = goi_gg_df,
    aes(x = RDA1, y = RDA2, label = goi_name),
    size = 4, col = "coral", max.overlaps = Inf
  ) +
  labs(col = "Genes of interest", fill = "SNP\ncount") +
  theme_void() 
dev.off()

pdf(here(plotpath, "goi_rda_by_chr.pdf"), width = 14, height = 5)
zscale = 8
ggplot() +
  geom_hex(data = rda_z, aes(x = RDA1, y = RDA2), binwidth = 0.2, alpha = 1) +
  stat_ellipse(data = rda_z, aes(x = RDA1, y = RDA2), level = 0.99, linetype = "dashed", size = 0.5) +
  geom_text(data = rda_varloadings, aes(x = RDA1 * zscale, y = RDA2 * zscale, label = var), size = 3, vjust = 1) +
  geom_segment(data = rda_varloadings, aes(x = 0, y = 0, xend = RDA1 * zscale, yend = RDA2 * zscale),
               arrow = arrow(length = unit(0.2, "cm"))) +
  scale_color_identity() +
  scale_x_continuous(expand = expansion(mult = 0.2)) +  # 5% margin
  scale_y_continuous(expand = expansion(mult = 0.2)) +
  scale_fill_gradientn(
    trans = "log10",
    colours = gray.colors(100, start = 1, end = 0)
  ) +
  geom_point(data = goi_gg_df, aes(x = RDA1, y = RDA2), col = "coral") +
  ggrepel::geom_text_repel(
    data = drop_na(gg_df, goi_name),
    aes(x = RDA1, y = RDA2, label = goi_name),
    size = 4, col = "coral", max.overlaps = Inf
  ) +
  labs(fill = "SNP count") +
  theme_void() +
  facet_wrap(~chr, scales = "free", nrow = 2) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), 
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt")
  )


ggplot() +
  #geom_hex(data = rda_sample, aes(x = RDA1, y = RDA2), bins = 50, alpha = 0.5) +
  geom_text(data = rda_varloadings, aes(x = RDA1 * zscale, y = RDA2 * zscale, label = var), size = 3, vjust = 1) +
  geom_segment(data = rda_varloadings, aes(x = 0, y = 0, xend = RDA1 * zscale, yend = RDA2 * zscale),
               arrow = arrow(length = unit(0.2, "cm"))) +
  scale_color_identity() +
  scale_fill_gradientn(
    trans = "log10",
    colours = gray.colors(100, start = 1, end = 0)
  ) +
  geom_point(data = goi_gg_df, aes(x = RDA1, y = RDA2), col = "coral") +
  ggrepel::geom_text_repel(
    data = drop_na(gg_df, goi_name),
    aes(x = RDA1, y = RDA2, label = goi_name),
    size = 4, col = "coral", max.overlaps = Inf
  ) +
  geom_point(data = top_rda, aes(x = RDA1, y = RDA2), col = "blue", size = 3) +
  ggrepel::geom_text_repel(
    data = top_rda,
    aes(x = RDA1, y = RDA2, label = locus),
    size = 5, col = "blue", max.overlaps = Inf
  ) +
  labs(fill = "SNP count") +
  theme_void() +
  facet_wrap(~chr, scales = "free", nrow = 2) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), 
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt")
  )
dev.off()
