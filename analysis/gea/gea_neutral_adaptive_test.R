# Read in heterozygosity for Adaptive SNPs

# Read in heterozygosity for Neutral SNPs

# OR SUBSET
# NEED TO THINK ABOUT WHETHER OR NOT TO LD PRUNE FIRST (probably should test both...)

# TEMP:
library(vcfR)
library(algatr)
library(tidyverse)
library(sf)
library(terra)
library(here)
source("general_functions.R")
vcf <- read.vcfR("data/processed_data/58-Sceloporus_JALMGF010000001.1.vcf.gz")
gen <- vcf[1:100000,]
coords <- 
  get_coords(sf = TRUE) %>% 
  filter(SampleID %in% colnames(vcf@gt[,-1]))
# CHECK THAT ORDER IS CORRECT
stopifnot(all(colnames(vcf@gt[,-1]) == coords$SampleID))

env <- rast(here("data", "env", "envstack.tif"))
env <- extract(env[[2:4]], coords, ID = FALSE)

dos <- vcf_to_dosage(gen)
dos_imputed <- simple_impute(dos)

mod <- rda_run(dos_imputed, env, coords, correctGEO = TRUE)

rda_sig_p <- rda_getoutliers(mod, naxes = "all", outlier_method = "p", p_adj = "fdr", sig = 0.01, plot = FALSE)
rda_snps <- rda_sig_p$rda_snps

dos_to_het <- function(x){
  het <- apply(x, 1, function(x) {
    sum(x == 1, na.rm = TRUE) / sum(!is.na(x))
  })
  return(het)
}

het_adaptive <- dos_to_het(dos_imputed[,rda_snps])

nrda <- length(rda_snps)

dos_neutral <- dos[,which(!(colnames(dos) %in% rda_snps))]
het_neutral <- map(1:1000, ~{
  dos_to_het(dos_neutral[,sample(1:ncol(dos_neutral), nrda)])
}) %>% 
list_transpose()

# For each individual, calculate a p-value
pvals <- map(names(het_adaptive), ~{
  value <- het_adaptive[.x]
  p_greater <- mean(het_neutral[[.x]] >= value)
  p_less <- mean(het_neutral[[.x]] <= value)
  p_value <- 2 * min(p_greater, p_less)
  return(c(greater = p_greater, less = p_less, both = p_value))
}) %>% bind_rows()

# 93% of the time, the individuals have significantly higher genetic diversity at adaptive SNPs than neutral SNPs
pvals %>%
  mutate_at(vars(everything()), ~p.adjust(.x, method = "fdr")) %>%
  mutate_at(vars(everything()), ~(.x < 0.05)) %>%
  summarise_all(mean)

# Calculate correlation between adaptive and neutral heterozygosity
t_het_neutral <- list_transpose(het_neutral)
rvals <- map_dbl(1:length(t_het_neutral), ~{
  cor(het_adaptive, t_het_neutral[[.x]], use = "pairwise.complete.obs")
}) 
# The mean correlation between adaptive and neutral heterozygosity is 0.67
mean(rvals)

mean_neutral <- map_dbl(het_neutral, ~mean(.x, na.rm = TRUE))
gg_df <-
  coords %>%
  mutate(neutral = mean_neutral, adaptive = het_adaptive) 

pdf(here("analysis", "gea", "het_adaptive_vs_neutral.pdf"), width = 5, height = 5)
ggplot(gg_df, aes(x = neutral, y = adaptive)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "black", lty = "dashed") +
  geom_smooth(method = "lm", se = TRUE) +
  labs(x = "Neutral SNP Heterozygosity", y = "Adaptive SNP Heterozygosity") +
  theme_classic()

gg_long <- pivot_longer(gg_df, cols = c(neutral, adaptive), names_to = "type", values_to = "het")
ggplot(gg_long) +
  geom_sf(aes(color = het)) +
  scale_color_viridis_c(option = "turbo") +
  theme_void() +
  facet_wrap(~type) +
  ggtitle("Heterozygosity at Adaptive and Neutral SNPs (Raw)")

gg_scaled <-
  gg_df %>%
  mutate(neutral = scale(neutral), adaptive = scale(adaptive)) %>%
  pivot_longer(cols = c(neutral, adaptive), names_to = "type", values_to = "het")
ggplot(gg_scaled) +
  geom_sf(aes(color = het)) +
  scale_color_viridis_c(option = "turbo") +
  theme_void() +
  facet_wrap(~type) +
  ggtitle("Heterozygosity at Adaptive and Neutral SNPs (Scaled)")

# Plot across space
dev.off()

# PROBLEM: does GEA make it such that SNPs with higher heterozygosity are favored because they are more likely to have correlations with the environment? Maybe we need to compare adaptive SNPs to neutral SNPs that have similar levels of heterozygosity at that SNP, but then compare heterozygosity for each individual (i.e., for SNPs with underlying levels of heterozygosity that are similar, do individuals have higher heterozygosity at adaptive SNPs than neutral SNPs?) Idk if this makes sense...
dos_neutral <- dos[,which(!(colnames(dos) %in% rda_snps))]
het_by_snp_neutral <- dos_to_het(t(dos_neutral))
het_by_snp_adaptive <- dos_to_het(t(dos_imputed[,rda_snps]))

names(rda_snps) <- rda_snps
het_neutral2 <- map(rda_snps, ~{
  heta <- round(het_by_snp_adaptive[.x], 2)
  hetn <- round(het_by_snp_neutral, 2) 
  idx <- which(hetn == heta)
  if (length(idx) < 100){
    warning("Not enough SNPs with similar heterozygosity")
    return(NULL)
  } else {
    idx <- names(idx)[1:100]
    return(dos_to_het(dos_neutral[,idx]))
  }
}) %>% compact() %>% list_transpose()

pvals2 <- map(names(het_adaptive), ~{
  value <- het_adaptive[.x]
  p_greater <- mean(het_neutral2[[.x]] >= value)
  p_less <- mean(het_neutral2[[.x]] <= value)
  p_value <- 2 * min(p_greater, p_less)
  return(c(greater = p_greater, less = p_less, both = p_value))
}) %>% bind_rows()

pvals2 %>%
  mutate_at(vars(everything()), ~p.adjust(.x, method = "fdr")) %>%
  mutate_at(vars(everything()), ~(.x < 0.05)) %>%
  summarise_all(mean)

# evaluate correaltions
t_het_neutral2 <- list_transpose(het_neutral2)
rvals2 <- map_dbl(1:length(t_het_neutral2), ~{
  cor(het_adaptive, t_het_neutral2[[.x]], use = "pairwise.complete.obs")
})
mean(rvals2)


mean_neutral <- map_dbl(het_neutral2, ~mean(.x, na.rm = TRUE))
gg_df <-
  coords %>%
  mutate(neutral = mean_neutral, adaptive = het_adaptive) 

pdf(here("analysis", "gea", "het_adaptive_vs_neutral2.pdf"), width = 5, height = 5)
ggplot(gg_df, aes(x = neutral, y = adaptive)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "black", lty = "dashed") +
  geom_smooth(method = "lm", se = TRUE) +
  labs(x = "Neutral SNP Heterozygosity", y = "Adaptive SNP Heterozygosity") +
  theme_classic()

gg_long <- pivot_longer(gg_df, cols = c(neutral, adaptive), names_to = "type", values_to = "het")
ggplot(gg_long) +
  geom_sf(aes(color = het)) +
  scale_color_viridis_c(option = "turbo") +
  theme_void() +
  facet_wrap(~type) +
  ggtitle("Heterozygosity at Adaptive and Neutral SNPs (Raw)")

gg_scaled <-
  gg_df %>%
  mutate(neutral = scale(neutral), adaptive = scale(adaptive)) %>%
  pivot_longer(cols = c(neutral, adaptive), names_to = "type", values_to = "het")
ggplot(gg_scaled) +
  geom_sf(aes(color = het)) +
  scale_color_viridis_c(option = "turbo") +
  theme_void() +
  facet_wrap(~type) +
  ggtitle("Heterozygosity at Adaptive and Neutral SNPs (Scaled)")

# Plot across space
dev.off()


# Compare with whole genome heterozygosity
source(here("analysis", "genetic_diversity", "genetic_diversity.R"))
het <- get_het()
gg_df <-
  coords %>% 
  left_join(het) %>%
  mutate(neutral = scale(Ho), adaptive = scale(het_adaptive))

pdf(here("analysis", "gea", "het_adaptive_vs_neutral3.pdf"), width = 5, height = 5)
ggplot(gg_df, aes(x = neutral, y = adaptive)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "black", lty = "dashed") +
  geom_smooth(method = "lm", se = TRUE) +
  labs(x = "Neutral SNP Heterozygosity (Scaled)", y = "Adaptive SNP Heterozygosity (Scaled)") +
  theme_classic()

gg_long <- pivot_longer(gg_df, cols = c(neutral, adaptive), names_to = "type", values_to = "het")
ggplot(gg_long) +
  geom_sf(aes(color = het)) +
  scale_color_viridis_c(option = "turbo") +
  theme_void() +
  facet_wrap(~type) +
  ggtitle("Heterozygosity at Adaptive and Neutral SNPs (Scaled)")
dev.off()

# IDEA:
1. Identify SNPs associated with adaptation to warmer temperatures
2. Predict where allele frequencies would need to increase to adapt to warmer temperatures
3. Compare the genetic diversity at these SNPs to the genetic diversity at neutral SNPs