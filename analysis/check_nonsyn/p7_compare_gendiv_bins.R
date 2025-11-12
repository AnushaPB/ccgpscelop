library(tidyverse)
library(here)
library(ggpubr)
library(sf)
library(terra)
plotpath <- here("analysis", "check_nonsyn", "plots")
source(here("general_functions.R"))

outpath <- here("analysis", "check_nonsyn", "outputs")

# 1. Get all het files
files_full <- list.files(outpath, pattern = "het_variants_MAF_.*\\.het$", full.names = TRUE)
files_nonsyn <- list.files(outpath, pattern = "het_nonsyn_variants_MAF_.*\\.het$", full.names = TRUE)
files_syn <- list.files(outpath, pattern = "het_syn_variants_MAF_.*\\.het$", full.names = TRUE)
all_files <- c(files_full, files_nonsyn, files_syn)

# 2. Read and combine
het_bins <- 
  all_files %>%
  set_names() %>%                    # keep filenames
  map_dfr(read_table, .id = "file") %>%  # read and stack
  mutate(
    # Label full vs nonsyn
    variant_set =
      case_when(
        str_detect(file, "het_variants_MAF_") ~ "all",
        str_detect(file, "het_nonsyn_variants_MAF_") ~ "non-synonymous",
        str_detect(file, "het_syn_variants_MAF_") ~ "synonymous",
        TRUE ~ NA_character_
      )
  ) %>%
  mutate(
    # Extract bin from filename: het_variants_MAF_0.00_0.01.het → 0.00_0.01
    bin = str_extract(file, "(?<=variants_MAF_)[0-9.]+_[0-9.]+"),
    # Optional: convert to pretty label 0.00–0.01
    bin = str_replace(bin, "_", "–")
  ) %>%
  mutate(Ho = 1 - `O(HOM)` / `N(NM)`) %>%  # Calculate observed heterozygosity
  select(bin, everything(), -file) %>%
  select(bin, SampleID = IID, Ho, variant_set, N = `N(NM)`)

het_bins %>% slice_head(n = 5)

# Calculate number of variants in each bin and variant set
variant_counts <-
  het_bins %>%
  group_by(bin, variant_set) %>%
  summarize(
    n_variants = mean(N, na.rm = TRUE)
  ) %>%
  ungroup()

head(variant_counts)

# Get model variables
model_df <- read_csv(here("analysis", "genetic_diversity", "outputs", "model_df.csv")) %>% rename(Ho_full = Ho)
syn_full <- read_table(here(outpath, "all_synonymous.het"), col_names = TRUE) %>% mutate(Ho_syn_full = 1 - (`O(HOM)` / `N(NM)`)) %>% select(SampleID = IID, Ho_syn_full)

het_bins_mod <- het_bins %>%
  left_join(model_df) %>%
  left_join(syn_full) %>%
  drop_na(Ho, Ho_full) %>%
  group_by(bin, variant_set) %>%
  group_modify(~{
    dat <- .x
    fit <- lm(Ho ~ Ho_full, data = dat)
    dat$Ho_resid <- resid(fit)
    fit <- lm(Ho ~ Ho_syn_full, data = dat)
    dat$Ho_resid_syn <- resid(fit)
    return(dat)
  }) %>%
  ungroup() %>%
  mutate(
    variant_set = factor(variant_set, levels = c("all", "synonymous", "non-synonymous"))
  )

  
alpha <- 0.05

library(broom)
# per-group tests
het_sig <- het_bins_mod %>%
  group_by(bin, variant_set) %>%
  group_modify(~{
    df <- .x

    res <- tibble(
      n = nrow(df),
      ho_tmean_slope = NA_real_,
      ho_tmean_p     = NA_real_,
      resid_tmean_slope = NA_real_,
      resid_tmean_p     = NA_real_
    )

    # Only run if there is variation in predictor
    if (nrow(df) >= 3 && dplyr::n_distinct(df$tmean_dif) >= 2) {
      # Ho ~ tmean_dif
      if (all(c("Ho","tmean_dif") %in% names(df)) && sum(!is.na(df$Ho)) >= 3) {
        fit1 <- try(lm(Ho ~ tmean_dif + csi_past + bio1 + Q, data = df), silent = TRUE)
        if (!inherits(fit1, "try-error")) {
          t1 <- broom::tidy(fit1)
          if ("tmean_dif" %in% t1$term) {
            res$ho_tmean_slope <- t1$estimate[t1$term == "tmean_dif"]
            res$ho_tmean_p     <- t1$p.value[t1$term == "tmean_dif"]
          }
        }
      }
      # Ho_resid ~ tmean_dif
      if (all(c("Ho_resid","tmean_dif") %in% names(df)) &&
          sum(!is.na(df$Ho_resid)) >= 3) {
        fit2 <- try(lm(Ho_resid ~ tmean_dif + csi_past + bio1 + Q, data = df), silent = TRUE)
        if (!inherits(fit2, "try-error")) {
          t2 <- broom::tidy(fit2)
          if ("tmean_dif" %in% t2$term) {
            res$resid_tmean_slope <- t2$estimate[t2$term == "tmean_dif"]
            res$resid_tmean_p     <- t2$p.value[t2$term == "tmean_dif"]
          }
        }
      }
    }
    res
  }) %>%
  ungroup() %>%
  # BH-adjust p-values within each variant_set across bins
  group_by(variant_set) %>%
  mutate(
    ho_tmean_p_adj    = ifelse(is.na(ho_tmean_p), NA_real_, p.adjust(ho_tmean_p, method = "fdr")),
    resid_tmean_p_adj = ifelse(is.na(resid_tmean_p), NA_real_, p.adjust(resid_tmean_p, method = "fdr")),
    ho_tmean_sig_raw    = !is.na(ho_tmean_p)    & ho_tmean_p    < alpha,
    resid_tmean_sig_raw = !is.na(resid_tmean_p) & resid_tmean_p < alpha,
    ho_tmean_sig_adj    = !is.na(ho_tmean_p_adj)    & ho_tmean_p_adj    < alpha,
    resid_tmean_sig_adj = !is.na(resid_tmean_p_adj) & resid_tmean_p_adj < alpha
  ) %>%
  ungroup()

het_bins_gg <- 
  het_bins_mod %>%
  drop_na(Ho) %>%
  left_join(het_sig, by = c("bin", "variant_set"))
  
het_bins_gg_filtered <-
  het_bins_gg %>%
  # keep only bins that exist within each variant_set
  semi_join(
    het_sig %>% filter(variant_set == "nonsyn") %>% select(bin),
    by = c("bin")
  )


alpha <- 0.05  # keep whatever you were using

pdf(here(plotpath, "Ho_by_MAF_bin.pdf"), width=15, height=6)

ggplot(het_bins_gg, aes(x = csi_past, y = Ho, col = bin, fill = bin)) +
  #geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  #stat_cor() +
  scale_color_viridis_d(option = "plasma") +
  scale_fill_viridis_d(option = "plasma") +
  scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "dotted")) +
  labs(x = make_pretty_names("csi_past"), y = "Binned Ho", color = "MAF Bin", fill = "MAF Bin") +
  facet_wrap(~variant_set, scales = "free_y")

ggplot(het_bins_gg, aes(x = tmean_dif, y = Ho, col = bin, fill = bin, linetype = ho_tmean_sig_adj)) +
  #geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  #stat_cor() +
  scale_color_viridis_d(option = "plasma") +
  scale_fill_viridis_d(option = "plasma") +
  #scale_alpha_manual(values = c("TRUE" = 0.8, "FALSE" = 0)) +
  scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "dotted")) +
  labs(x = make_pretty_names("tmean_dif"), y = "Binned Ho", color = "MAF Bin", fill = "MAF Bin") +
  facet_wrap(~variant_set, scales = "free_y")

ggplot(het_bins_gg, aes(x = tmean_dif, y = Ho_resid, col = bin, fill = bin, linetype = resid_tmean_sig_adj)) +
  #geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  #stat_cor() +
  scale_color_viridis_d(option = "plasma") +
  scale_fill_viridis_d(option = "plasma") +
  #scale_alpha_manual(values = c("TRUE" = 0.8, "FALSE" = 0)) +
  scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "dotted")) +
  labs(x = make_pretty_names("tmean_dif"), y = "Residual(binnned ~ whole genome Ho)", color = "MAF Bin", fill = "MAF Bin") +
  facet_wrap(~variant_set, scales = "free_y")

dev.off()

gls_data <-
  het_bins_mod %>%
  group_by(x, y, variant_set, bin) %>%
  summarize(
    across(c(Ho_resid, Ho_resid_syn, Ho, tmean_dif, csi_past, bio1, Q), ~ mean(.x, na.rm = TRUE)
  )) %>%
  st_as_sf(coords = c("x", "y"), crs = 4326) %>%
  st_transform(3310) %>%
  bind_cols(st_coordinates(.)) %>%
  drop_na() %>%
  group_by(bin, variant_set) %>%
  group_split()

library(nlme)
run_gls <- function(data, y) {
  f <- formula(paste0(y, " ~ tmean_dif + csi_past + bio1 + Q"))
  gls_model <- gls(f, correlation = corExp(form = ~ X + Y, nugget = FALSE), data = data) 
  gls_tb <- 
    summary(gls_model)$tTable %>% 
    data.frame() %>% 
    rownames_to_column(var = "term") %>%
    mutate(bin = unique(data$bin), 
    variant_set = unique(data$variant_set)) 
  
  fitted_tb <- 
    get_fitted_tmean_values(gls_model, data) %>%
    tibble(
      bin       = unique(data$bin),
      variant_set = unique(data$variant_set)
    )

  return(list(gls = gls_tb, fitted = fitted_tb))
}

get_fitted_tmean_values <- function(gls_model, data) {
  xseq <- seq(min(data$tmean_dif, na.rm = TRUE),
              max(data$tmean_dif, na.rm = TRUE),
              length.out = 100)

  nd <- tibble(
    tmean_dif = xseq,
    csi_past  = mean(data$csi_past, na.rm = TRUE),
    bio1      = mean(data$bio1,   na.rm = TRUE),
    Q         = mean(data$Q,      na.rm = TRUE),
    # coordinates are only for correlation structure; set to group means
    X = mean(data$X, na.rm = TRUE),
    Y = mean(data$Y, na.rm = TRUE)
  )

  pred <- predict(gls_model, newdata = nd)

  tb <-
    tibble(
      tmean_dif = xseq,
      fitted    = pred
    )

  return(tb)
}

gls_result_Ho <- map(gls_data, ~run_gls(.x, "Ho"), .progress = TRUE)  
gls_table_Ho <- list_transpose(gls_result_Ho)[["gls"]] %>% bind_rows() %>% mutate(response = "Ho")
gls_fitted_Ho <- list_transpose(gls_result_Ho)[["fitted"]] %>% bind_rows() %>% mutate(response = "Ho")

gls_result_Ho_resid <- map(gls_data, ~run_gls(.x, "Ho_resid"), .progress = TRUE) 
gls_table_Ho_resid <- list_transpose(gls_result_Ho_resid)[["gls"]] %>% bind_rows() %>% mutate(response = "Ho_resid")
gls_fitted_Ho_resid <- list_transpose(gls_result_Ho_resid)[["fitted"]] %>% bind_rows() %>% mutate(response = "Ho_resid")

gls_result_Ho_resid_syn <- map(gls_data, ~run_gls(.x, "Ho_resid_syn"), .progress = TRUE)
gls_table_Ho_resid_syn <- list_transpose(gls_result_Ho_resid_syn)[["gls"]] %>% bind_rows() %>% mutate(response = "Ho_resid_syn")
gls_fitted_Ho_resid_syn <- list_transpose(gls_result_Ho_resid_syn)[["fitted"]] %>% bind_rows() %>% mutate(response = "Ho_resid_syn")

variant_counts %>%
  arrange(n_variants)

gls_table_result <-
  bind_rows(gls_table_Ho_resid, gls_table_Ho_resid_syn, gls_table_Ho) %>%
  mutate(p_adj = p.adjust(p.value, method = "fdr")) %>%
  mutate(significant = ifelse(p_adj < alpha, TRUE, FALSE)) %>%
  left_join(variant_counts) %>%
  mutate(
    variant_set = factor(variant_set, levels = c("all", "synonymous", "non-synonymous"))
  ) 

gls_fitted_result <-
  bind_rows(gls_fitted_Ho_resid, gls_fitted_Ho)

n_variant_under_1k <- 
  variant_counts %>% 
  filter(n_variants < 1000) %>%
  arrange(bin) %>%
  head(1) %>% 
  pull(bin)

bins_with_over_1k <- 
  variant_counts %>%
  filter(variant_set == "non-synonymous") %>%
  filter(n_variants >= 1000) %>%
  arrange(bin) %>%
  pull(bin)


color_values <- c("non-synonymous" = "red2", "synonymous" = "blue2", "all" = "#5f5f5f")

pdf(here(plotpath, "gls_Ho_by_MAF_bin.pdf"), width=10, height=5)
ggplot(filter(gls_table_result, response == "Ho", term == "tmean_dif")) +
  geom_rect(aes(xmin = n_variant_under_1k, xmax = Inf, ymin = -Inf, ymax = Inf), fill = "#dbe7ff", inherit.aes = FALSE, alpha = 0.5) +
  geom_text(aes(x = n_variant_under_1k, y = min(Value, na.rm = TRUE), label = " < 1000 non-synonymous variants", hjust = 0, vjust = -1), color = "#0b0039", inherit.aes = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_point(aes(x = bin, y = Value, shape = significant, color = variant_set, size = (n_variants))) +
  geom_line(aes(x = bin, y = Value, group = variant_set, color = variant_set)) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1)) +
  scale_color_manual(values = color_values) +
  scale_size_continuous(trans = "log10") + 
  theme_classic() +
  # rotate axis labels
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = "MAF Bin", y = "GLS coefficient", 
    color = "Variant Set", shape = "Significant", size = "Number of Variants",
    title = "Effect of recent climate change on Ho"
  )

ggplot(filter(gls_table_result, response == "Ho_resid", term == "tmean_dif")) +
  geom_rect(aes(xmin = n_variant_under_1k, xmax = Inf, ymin = -Inf, ymax = Inf), fill = "#dbe7ff", inherit.aes = FALSE, alpha = 0.5) +
  geom_text(aes(x = n_variant_under_1k, y = min(Value, na.rm = TRUE), label = " < 1000 non-synonymous variants", hjust = 0, vjust = -1), color = "#0b0039", inherit.aes = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_point(aes(x = bin, y = Value, shape = significant, color = variant_set, size = (n_variants))) +
  geom_line(aes(x = bin, y = Value, group = variant_set, color = variant_set)) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1)) +
  scale_color_manual(values = color_values) +
  scale_size_continuous(trans = "log10") + 
  theme_classic() +
  # rotate axis labels
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = "MAF Bin", y = "GLS coefficient", 
    color = "Variant Set", shape = "Significant", size = "Number of Variants",
    title = "Effect of recent climate change on residual(binned Ho ~ genome-wide Ho)"
  ) 

ggplot(filter(gls_table_result, response == "Ho_resid_syn", term == "tmean_dif", variant_set != "all")) +
  geom_rect(aes(xmin = n_variant_under_1k, xmax = Inf, ymin = -Inf, ymax = Inf), fill = "#dbe7ff", inherit.aes = FALSE, alpha = 0.5) +
  geom_text(aes(x = n_variant_under_1k, y = min(Value, na.rm = TRUE), label = " < 1000 non-synonymous variants", hjust = 0, vjust = -1), color = "#0b0039", inherit.aes = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_point(aes(x = bin, y = Value, shape = significant, color = variant_set, size = (n_variants))) +
  geom_line(aes(x = bin, y = Value, group = variant_set, color = variant_set)) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1)) +
  scale_color_manual(values = color_values) +
  scale_size_continuous(trans = "log10") + 
  theme_classic() +
  # rotate axis labels
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = "MAF Bin", y = "GLS coefficient", 
    color = "Variant Set", shape = "Significant", size = "Number of Variants",
    title = "Effect of recent climate change on residual(binned Ho ~ genome-wide synonymous Ho)"
  ) 
dev.off()

# TESTING RESIDUALS WITHIN BINS ------------------------------------------
gls_data2 <-
  gls_data %>%
  bind_rows() %>%
  select(-Ho_resid, -Ho_resid_syn) %>%
  pivot_wider(names_from = variant_set, values_from = Ho) %>%
  mutate(
    Ho_resid_bin = resid(lm(`non-synonymous` ~ all)),
    Ho_resid_syn_bin = resid(lm(`non-synonymous` ~ synonymous))
  ) %>% 
  group_by(bin) %>%
  group_split()

gls_result_Ho_resid_syn_bin <- map(gls_data2, ~run_gls(.x, "Ho_resid_syn_bin"), .progress = TRUE)
gls_table_Ho_resid_syn_bin <- list_transpose(gls_result_Ho_resid_syn_bin)[["gls"]] %>% bind_rows() %>% mutate(response = "Ho_resid_syn_bin")
gls_result_Ho_resid_bin <- map(gls_data2, ~run_gls(.x, "Ho_resid_bin"), .progress = TRUE)
gls_table_Ho_resid_bin <- list_transpose(gls_result_Ho_resid_bin)[["gls"]] %>% bind_rows() %>% mutate(response = "Ho_resid_bin")

gls_table_Ho_bin <- 
  bind_rows(gls_table_Ho_resid_bin, gls_table_Ho_resid_syn_bin) %>%
  mutate(p_adj = p.adjust(p.value, method = "fdr")) %>%
  mutate(significant = ifelse(p_adj < alpha, TRUE, FALSE)) %>%
  left_join(variant_counts %>% filter(variant_set == "non-synonymous")) %>%
  filter(term %in% c("bio1", "csi_past", "tmean_dif"))

pdf(here(plotpath, "gls_Ho_resid_by_MAF_bin.pdf"), width=10, height=5)
ggplot(filter(gls_table_Ho_bin, response == "Ho_resid_syn_bin", term == "tmean_dif")) +
  geom_rect(aes(xmin = n_variant_under_1k, xmax = Inf, ymin = -Inf, ymax = Inf), fill = "#dbe7ff", inherit.aes = FALSE, alpha = 0.5) +
  geom_text(aes(x = n_variant_under_1k, y = min(Value, na.rm = TRUE), label = " < 1000 non-synonymous variants", hjust = 0, vjust = -1), color = "#0b0039", inherit.aes = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_point(aes(x = bin, y = Value, shape = significant, size = (n_variants), color = variant_set)) +
  geom_line(aes(x = bin, y = Value, group = 1, color = variant_set)) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1)) +
  scale_color_manual(values = color_values) +
  scale_size_continuous(trans = "log10") + 
  theme_classic() +
  # rotate axis labels
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = "MAF Bin", y = "GLS coefficient", 
    shape = "Significant", size = "Number of Variants", color = "Variant Set", 
    title = "Effect of recent climate change on residual(binned non-synonymous Ho ~ binned synonymous Ho)"
  ) 


ggplot(filter(gls_table_Ho_bin, response == "Ho_resid_syn_bin")) +
  geom_rect(aes(xmin = n_variant_under_1k, xmax = Inf, ymin = -Inf, ymax = Inf), fill = "#dbe7ff", inherit.aes = FALSE, alpha = 0.5) +
  geom_text(aes(x = n_variant_under_1k, y = min(Value, na.rm = TRUE), label = " < 1000 non-synonymous variants", hjust = 0, vjust = -1), color = "#0b0039", inherit.aes = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_point(aes(x = bin, y = Value, shape = significant, size = (n_variants), color = make_pretty_names(term))) +
  geom_line(aes(x = bin, y = Value, group = term, color = make_pretty_names(term))) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1)) +
  scale_size_continuous(trans = "log10") + 
  theme_classic() +
  # rotate axis labels
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = "MAF Bin", y = "GLS coefficient", 
    shape = "Significant", size = "Variable", color = "Variant Set", 
    title = "Effect on residual(binned non-synonymous Ho ~ binned synonymous Ho)"
  ) 

dev.off()

# MAPPING -------------------------------------------------------
het_bins_sf <-
  het_bins_mod %>%
  st_as_sf(coords = c("x", "y"), crs = 4326) %>%
  st_transform(3310) 

bins <- unique(het_bins_sf %>% pull(bin))
# Select evenly 5 bins
subbins <- bins[seq(1, length(bins), length.out = 6)]


mapping_data <-
  gls_data2 %>% 
  bind_rows()  %>%
  filter(bin %in% subbins)

ca <- get_ca()
plt1 <-
  ggplot() +
  geom_sf(data =ca) +
  geom_sf(data = mapping_data, aes(color = all), size = 2) +
  scale_color_viridis_c(option = "plasma") +
  theme_void() +
  facet_wrap(~bin, nrow = 1) +
  theme(strip.text = element_text(size = 14), title = element_text(size = 18)) + 
  ggtitle("All Variants Ho")

plt2 <-
  ggplot() +
  geom_sf(data =ca) +
  geom_sf(data = mapping_data, aes(color = `non-synonymous`), size = 2) +
  scale_color_viridis_c(option = "plasma") +
  theme_void() +
  facet_wrap(~bin, nrow = 1) +
  theme(strip.text = element_text(size = 14), title = element_text(size = 18)) +
  ggtitle("Nonsynonymous Variants Ho")

plt3 <-
  ggplot() +
  geom_sf(data =ca) +
  geom_sf(data = mapping_data, aes(color = synonymous), size = 2) +
  scale_color_viridis_c(option = "plasma") +
  theme_void() +
  facet_wrap(~bin, nrow = 1) +
  theme(strip.text = element_text(size = 14), title = element_text(size = 18)) +
  ggtitle("Synonymous Variants Ho")

plt4 <-
  ggplot() +
  geom_sf(data =ca) +
  geom_sf(data = mapping_data, aes(color = Ho_resid_syn_bin), size = 2) +
  scale_color_viridis_c(option = "plasma") +
  theme_void() +
  facet_wrap(~bin, nrow = 1) +
  theme(strip.text = element_text(size = 14), title = element_text(size = 18)) +
  ggtitle("Binned residual(non-synonymous ~ synonymous)")

library(cowplot)
png(here(plotpath, "map_Ho_binned.png"), width=4 * length(subbins)*300, height=5*4*300, , res = 300)
plot_grid(plt1, plt2, plt3, plt4, ncol = 1)
dev.off()

model_sf <- 
  model_df %>%
  st_as_sf(coords = c("x", "y"), crs = 4326) %>%
  st_transform(3310)

png(here(plotpath, "genome_wide_Ho_map.png"), width=6*300, height=6*300, , res = 300)
ggplot() +
  geom_sf(data =ca) +
  geom_sf(data = model_sf, aes(color = Ho_full), size = 2) +
  scale_color_viridis_c(option = "plasma") +
  theme_void() +
  labs(col = "Genome-wide Ho")
dev.off()


# GIF
library(gganimate)
library(ggplot2)
library(sf)
library(viridis)


all_variants_maf <- 
  het_bins_sf %>%
  filter(variant_set == "all") 

ca <- get_ca()
p <-
  ggplot() +
  geom_sf(data = ca) +
  geom_sf(data = all_variants_maf, aes(color = Ho), size = 2, alpha = 1) +
  scale_color_viridis_c(option = "plasma") +
  theme_void() +
  labs(title = "All Variants Ho — MAF Bin: {current_frame}") +
  transition_manual(bin)  +
  theme(legend.position = "bottom")

anim <- animate(
  p,
  fps = 10,                  # <-- increase fps to make it faster
  width = 5*300,
  height = 5*300,
  res = 300,   
  renderer = gifski_renderer()
)

# Save as gif
anim_save(here(plotpath, "all_variants_Ho.gif"), animation = anim)


# PLOTTING FITTED VALUES -------------------------------------------------------

sig_bin <- 
  gls_table_result %>% 
  filter(term == "tmean_dif") %>%
  select(significant, bin, variant_set, response) %>%
  mutate(response = paste0(response, "_significant")) %>%
  pivot_wider(names_from = response, values_from = significant) 

gls_gg <-
  gls_fitted_result %>%
  left_join(sig_bin) %>%
  mutate(
    variant_set = factor(variant_set, levels = c("all", "synonymous", "non-synonymous"))
  )  
names(gls_gg)

pdf(here(plotpath, "GLS_Ho_by_MAF_bin_fitted.pdf"), width=15, height=6)
ggplot(gls_gg, aes(x = tmean_dif, y = fitted, col = bin, fill = bin, linetype = Ho_significant)) +
  #geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  #stat_cor() +
  scale_color_viridis_d(option = "plasma") +
  scale_fill_viridis_d(option = "plasma") +
  #scale_alpha_manual(values = c("TRUE" = 0.8, "FALSE" = 0)) +
  scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "dotted")) +
  labs(x = make_pretty_names("tmean_dif"), y = "Binned Ho", color = "MAF Bin", fill = "MAF Bin") +
  facet_wrap(~variant_set, scales = "free_y")

ggplot(gls_gg, aes(x = tmean_dif, y = fitted, col = bin, fill = bin, linetype = Ho_resid_significant)) +
  #geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  #stat_cor() +
  scale_color_viridis_d(option = "plasma") +
  scale_fill_viridis_d(option = "plasma") +
  #scale_alpha_manual(values = c("TRUE" = 0.8, "FALSE" = 0)) +
  scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "dotted")) +
  labs(x = make_pretty_names("tmean_dif"), y = "Residual(binnned ~ whole genome Ho)", color = "MAF Bin", fill = "MAF Bin") +
  facet_wrap(~variant_set, scales = "free_y")

ggplot(filter(gls_gg, Ho_significant), aes(x = tmean_dif, y = fitted, col = bin, fill = bin, linetype = Ho_significant)) +
  #geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  stat_cor() +
  scale_color_viridis_d(option = "plasma") +
  scale_fill_viridis_d(option = "plasma") +
  #scale_alpha_manual(values = c("TRUE" = 0.8, "FALSE" = 0)) +
  scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "dotted")) +
  labs(x = make_pretty_names("tmean_dif"), y = "Binned Ho", color = "MAF Bin", fill = "MAF Bin") +
  facet_wrap(~variant_set, scales = "free_y")

ggplot(filter(gls_gg, Ho_resid_significant), aes(x = tmean_dif, y = fitted, col = bin, fill = bin, linetype = Ho_resid_significant)) +
  #geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  stat_cor() +
  scale_color_viridis_d(option = "plasma") +
  scale_fill_viridis_d(option = "plasma") +
  #scale_alpha_manual(values = c("TRUE" = 0.8, "FALSE" = 0)) +
  scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "dotted")) +
  labs(x = make_pretty_names("tmean_dif"), y = "Residual(binnned ~ whole genome Ho)", color = "MAF Bin", fill = "MAF Bin") +
  facet_wrap(~variant_set, scales = "free_y")
dev.off()