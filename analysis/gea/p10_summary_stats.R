library(tidyverse)
library(here)
source(here("general_functions.R"))
source(here("analysis", "gea", "functions_selection_stats.R"))
outpath <- here("analysis", "gea", "outputs")

# Get non-synonymous GEA variants in genes 
genes_nonsyn <- read_csv(here(outpath, "bio1ndvi_gea_genes_nonsyn.csv"))

# Get all annotated genes
all_genes <-read_csv(here("analysis", "gea", "outputs", "all_genes_list.csv"))

# Filter to only include genes with non-synonymous GEA variants
genes_org <- 
  all_genes %>% 
  # TODO: this shouldn't reall need to be two steps
  filter(full_name %in% genes_nonsyn$full_name) %>%
  drop_na(uniprot_id)

print(paste("Number of GEA SNPs in genes:", nrow(genes_nonsyn))) # 17,731
print(paste("Number of unique genes with GEA SNPs:", length(unique(genes_nonsyn$full_name)))) # 11,671
print(paste("Number of unique GEA UniprotIDs:", nrow(genes_org %>% distinct(uniprot_id)))) # 3,442
print(paste("Number of all genes with UniprotIDs:", nrow(all_genes))) # 30,349
print(paste("Number of all unique UniprotIDs:", nrow(all_genes %>% distinct(uniprot_id)))) # 18,357


# PLOTTING RDA
# IMPORTANT: NEED TO DO BY SCAFFOLD!!!!!!
rda_z <- read_csv(here("analysis", "gea", "outputs", "58-Sceloporus_unscaledloadings.csv"))
# Filter to just locus scores
rda_z <- rda_z %>% filter(score == "species") %>% rename(locus = label)
rda_varloadings <- 
  read_csv(here("analysis", "gea", "outputs", "58-Sceloporus_biplot.csv")) %>% 
  mutate(var = make_pretty_names(var)) %>%
  mutate(var = gsub("Contemporary t", "T", var))

# Filter to just chromosomes 
rda_z_chr <- 
  rda_z %>% 
  filter(grepl("chr", scaff)) %>%
  #filter(scaff %in% paste0("chr", 1:5)) %>%
  mutate(scaff = factor(scaff, levels = paste0("chr", 1:11))) %>%
  # RESCALING FOR PLOTTING
  group_by(scaff) %>%
  mutate(across(starts_with("RDA"),  ~ as.vector(scale(.))))

rda_varloadings <- 
  rda_varloadings %>% 
  filter(grepl("chr", scaff))%>%
  #filter(scaff %in% paste0("chr", 1:5)) %>%
  mutate(scaff = factor(scaff, levels = paste0("chr", 1:11))) 

# Get RDA loadings for significant loci
rda_sig_z <- 
  # DONT USE LOADINGS FROM HERE (they are from adaptive rda)
  read_csv(here("analysis", "adaptive", "outputs", "RDA_unscaledloadings_58-Sceloporus_bio1ndvi_gea_mod.csv")) %>%
  dplyr::select(locus = label) %>%
  distinct(locus) %>%
  left_join(rda_z_chr, by = "locus") 

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
  rda_z_chr %>%
  group_by(scaff) %>%
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
  rda_z_chr %>% 
  group_by(scaff) %>% 
  slice_max(RDA1, n = 100, with_ties = FALSE) %>% 
  ungroup()
top_rda2 <- 
  rda_z_chr %>% 
  group_by(scaff) %>% 
  slice_max(RDA2, n = 100, with_ties = FALSE) %>% 
  ungroup()
top_rda <- bind_rows(top_rda1, top_rda2) %>% distinct()
genes_nonsyn %>% filter(locus %in% top_rda$locus) %>% drop_na(gene_name) %>% dplyr::select(gene_name) %>% print(n = nrow(.)) # 1


# Top 5 for plotting
top_rda1 <- 
  rda_z_chr %>% 
  group_by(scaff) %>% 
  slice_max(RDA1, n = 5, with_ties = FALSE) %>% 
  ungroup()
top_rda2 <- 
  rda_z_chr %>% 
  group_by(scaff) %>% 
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

pdf(here(plotpath, "goi_rda.pdf"), width = 6, height = 5)
zscale = 8

ggplot() +
  #geom_point(data = rda_sample, aes(x = RDA1, y = RDA2), col = "gray", alpha = 0.5) +
  geom_hex(data = rda_z_chr, aes(x = RDA1, y = RDA2), binwidth = 0.2, alpha = 1) +
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
  geom_hex(data = rda_z_chr, aes(x = RDA1, y = RDA2), binwidth = 0.2, alpha = 1) +
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
  geom_hex(data = rda_z_chr, aes(x = RDA1, y = RDA2), binwidth = 0.2, alpha = 1) +
  stat_ellipse(data = rda_z_chr, aes(x = RDA1, y = RDA2), level = 0.99, linetype = "dashed", size = 0.5) +
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
  facet_wrap(~scaff, scales = "free", nrow = 2) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), 
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt")
  )


ggplot() +
  geom_hex(data = rda_sample, aes(x = RDA1, y = RDA2), bins = 50, alpha = 0.5) +
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
  facet_wrap(~scaff, scales = "free", nrow = 2) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), 
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt")
  )
dev.off()
