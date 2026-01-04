library(tidyverse)
library(here)
library(emmeans)
source(here("general_functions.R"))
outpath <- here("analysis", "melanization", "outputs")

full_pc <- read_tsv("outputs/pca.eigenvec") %>% select(SampleID = IID, genetic_PC1 = PC1, genetic_PC2 = PC2)

csq <- read_csv(here(outpath, "csq.csv"))
nonsyn <- csq %>% filter(grepl("missense|stop_gained|start_lost", csq))

dos <- read_table(here(outpath, "melanization_genes.raw"))
dos_df <- dos %>% rename(SampleID = IID) %>% select(starts_with("chr"), SampleID)
snp_names <- dos_df %>% select(starts_with("chr")) %>% colnames()
names(snp_names) <- snp_names

# parse scaffold + pos from your plink-style SNP IDs
nonsyn_snps <- 
  data.frame(snp = snp_names) %>%
  mutate(
    scaffold = str_extract(snp, "^[^_]+"),                       # "chr2"
    position = as.numeric(str_extract(snp, "(?<=_)[0-9]+"))      # 97054864
  ) %>% 
  inner_join(nonsyn)

env_df <- 
  read_csv(here("analysis", "genetic_diversity", "outputs", "model_df.csv")) %>%
  select(SampleID, bio1, NDVI, fire_frq, evt)

mod_df <- 
  left_join(env_df, dos_df, by = "SampleID") %>%
  mutate(fire_frq = factor(fire_frq, levels = c("Infrequent/No Fire Regime", "Intermediate Frequency", "High Frequency"))) %>%
  left_join(full_pc)

# Set fire_frq reference level to Infrequent/No Fire Regime
mod_df$fire_frq <- relevel(mod_df$fire_frq, ref = "Infrequent/No Fire Regime")

nonsyn_snpnames <- nonsyn_snps$snp
names(nonsyn_snpnames) <- nonsyn_snpnames
models <- 
  map(nonsyn_snpnames, ~{
    f <- paste0(.x, " ~ bio1 + NDVI + fire_frq + evt + genetic_PC1 + genetic_PC2")
    mod <- lm(f, data = mod_df)
  }, .progress = TRUE)

fire_results <-
  map(models, ~{
    emmeans(.x, ~ fire_frq) %>% pairs(infer = TRUE) %>% data.frame()
  }, .progress = TRUE) 

fire_results_df <- bind_rows(fire_results, .id = "snp")

fire_sig <- 
  fire_results_df %>% 
  filter(p.value < 0.05)

fire_sig

fire_sig %>% arrange(desc(abs(estimate))) %>% head()

source(here("analysis", "gea", "functions_selection_stats.R"))
melanization_locations <- 
  read_tsv(here("analysis", "melanization", "outputs", "melanization_genes.bed"), col_names = c("scaffold", "start", "end", "gene")) %>%
  mutate(start = start + 1) %>% # Convert from bed to gff
  left_join(get_all_genes_bed())
library(tidyverse)

# make sure intervals are numeric
melanization_locations2 <- melanization_locations %>%
  mutate(
    start = as.numeric(start),
    end   = as.numeric(end)
  )

# parse scaffold + pos from your plink-style SNP IDs
fire_sig_parsed <- fire_sig %>%
  mutate(
    scaffold = str_extract(snp, "^[^_]+"),                       # "chr2"
    pos      = as.numeric(str_extract(snp, "(?<=_)[0-9]+"))      # 97054864
  )

# annotate each SNP with the matching gene record (full_name)
fire_sig_annotated <- fire_sig_parsed %>%
  mutate(
    full_name = purrr::map2_chr(scaffold, pos, \(scaf, p) {
      if (is.na(scaf) || is.na(p)) return(NA_character_)
      hit <- melanization_locations2 %>%
        filter(scaffold == scaf, start <= p, end >= p) %>%
        pull(full_name)
      if (length(hit) == 0) NA_character_ else hit[[1]]
    }),
    gene = str_extract(full_name, "gene=([^;]+)") %>% str_remove("^gene=")
  )

fire_sig_annotated %>% arrange(desc(abs(estimate))) %>% select(snp, gene, estimate)
fire_sig_annotated %>% group_by(gene) %>% count()

mc1r <- fire_sig_annotated %>% filter(gene == "MC1R") 
mc1r_df <- mod_df %>% select(all_of(mc1r$snp), SampleID, fire_frq) 
mc1r_pca <- prcomp(mc1r_df[, mc1r$snp], center = TRUE, scale. = TRUE)
mc1r_df$PC1 <- mc1r_pca$x[, 1]
mc1r_df$PC2 <- mc1r_pca$x[, 2]
library(sf)
mc1r_df <- left_join(get_coords(sf = TRUE), mc1r_df, by = "SampleID")
mc1r_df_long <- mc1r_df %>% st_drop_geometry() %>% pivot_longer(starts_with("chr"))
pdf(here("analysis", "melanization", "plots", "mc1r_pca.pdf"))
ggplot(mc1r_df) +
  geom_sf(aes(col = PC1)) +
  scale_color_viridis_c(option = "turbo")

ggplot(mc1r_df) +
  geom_sf(aes(col = fire_frq))+
  scale_color_viridis_d(option = "turbo")

ggplot(mc1r_df) +
  geom_boxplot(aes(x = fire_frq, y = PC1))
dev.off()
# RDA --------------------------------------

library(tidyverse)
library(here)
library(vegan)

# --- load your objects (from your snippet) ---
outpath <- here("analysis", "melanization", "outputs")

dos <- read_table(here(outpath, "melanization_genes.raw"))
dos_df <- dos %>% rename(SampleID = IID) %>% select(starts_with("chr"), SampleID)

env_df <-
  read_csv(here("analysis", "genetic_diversity", "outputs", "model_df.csv")) %>%
  select(SampleID, bio1, NDVI, fire_frq, evt)

mod_df <-
  left_join(env_df, dos_df, by = "SampleID") %>%
  mutate(fire_frq = factor(fire_frq,
                           levels = c("Infrequent/No Fire Regime",
                                      "Intermediate Frequency",
                                      "High Frequency")))

mod_df$fire_frq <- relevel(mod_df$fire_frq, ref = "Infrequent/No Fire Regime")

# --- build Y (dosage matrix) ---
snp_names <- dos_df %>% select(starts_with("chr")) %>% names()

Y <- 
  mod_df %>%
  select(all_of(snp_names)) %>%
  mutate(across(everything(), as.numeric)) %>%
  as.matrix()

# optional: handle missing genotypes (simple mean impute per SNP)
library(algatr)
Y <- simple_impute(Y)

# Drop columns with all NA values
Y <- Y[, colSums(is.na(Y)) < nrow(Y)] # CHECK THIS

# recommended: scale SNPs so high-variance SNPs don't dominate
Y <- scale(Y, center = TRUE, scale = TRUE)

# --- build predictor dataframe aligned to Y rows ---
dfX <- mod_df %>%
  select(SampleID, fire_frq, bio1, NDVI, evt) %>%
  drop_na()

# keep only complete-case rows for predictors and match Y accordingly
keep <- match(dfX$SampleID, mod_df$SampleID)
Y <- Y[keep, , drop = FALSE]

# --- run RDA (categorical predictor fire_frq included) ---
library(vegan)
fit <- rda(Y ~ fire_frq + bio1 + NDVI + evt, data = dfX)

# # tests
# anova_overall <- anova.cca(fit, permutations = 999)
# anova_terms   <- anova.cca(fit, by = "term", permutations = 999)
# anova_axes    <- anova.cca(fit, by = "axis", permutations = 999)

# print(anova_overall)
# print(anova_terms)
# print(anova_axes)

# # quick summary + site scores if you want to plot
# sum_fit <- summary(fit)
# site_scores <- as.data.frame(scores(fit, display = "sites"))
# site_scores$SampleID <- dfX$SampleID
# site_scores$fire_frq <- dfX$fire_frq


