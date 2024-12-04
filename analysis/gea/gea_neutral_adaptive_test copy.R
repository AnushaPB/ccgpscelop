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
source(here("general_functions.R"))
vcf <- read.vcfR(here("data/ccgp_data/58-Sceloporus_pruned_mil.vcf.gz"))

gen <- vcf[1:100000,]
dos <- vcf_to_dosage(gen)
dos_imputed <- simple_impute(dos)
rm(vcf)

coords <- 
  get_coords(sf = TRUE) %>% 
  filter(SampleID %in% colnames(vcf@gt[,-1]))
# CHECK THAT ORDER IS CORRECT
stopifnot(all(colnames(vcf@gt[,-1]) == coords$SampleID))

env <- rast(here("data", "env", "envstack.tif"))[[2:4]]
env <- rast(here("data", "env", "california_env_stack.tif"))
env <- env[[c(1,5)]]
names(env) <- c("BIO1", "BIO5")
env <- extract(env, coords, ID = FALSE)

mod <- rda_run(dos_imputed, env, coords, correctGEO = TRUE)

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

# PROBLEM: does GEA make it such that SNPs with higher heterozygosity are favored because they are more likely to have correlations with the environment? 
# IDEA: do a maf filter with the cutoff being the minimum maf of the adaptive SNPs. Idk if this is biasing in the opposite direction though...
freq_adaptive <- apply(dos[,rda_snps], 2, function(x) mean(x, na.rm = TRUE)/2)
freq_min <- min(freq_adaptive)
freq_max <- max(freq_adaptive)

dos_neutral <- dos[,which(!(colnames(dos) %in% rda_snps))]
freq_neutral <- apply(dos_neutral, 2, function(x) mean(x, na.rm = TRUE)/2)
dos_neutral_maf <- dos_neutral[, freq_neutral >= freq_min & freq_neutral <= freq_max]
# Check that the MAFs are the same
freq_neutral_maf <- apply(dos_neutral_maf, 2, function(x) mean(x, na.rm = TRUE)/2)
stopifnot(all(range(freq_neutral_maf) == range(freq_adaptive)))

summary(freq_adaptive)
summary(freq_neutral_maf)

het_neutral2 <- map(1:1000, ~{
  dos_to_het(dos_neutral_maf[,sample(1:ncol(dos_neutral_maf), nrda)])
}) %>% 
list_transpose()

# For each individual, calculate a p-value
pvals <- map(names(het_adaptive), ~{
  value <- het_adaptive[.x]
  p_greater <- mean(het_neutral2[[.x]] >= value)
  p_less <- mean(het_neutral2[[.x]] <= value)
  p_value <- 2 * min(p_greater, p_less)
  return(c(greater = p_greater, less = p_less, both = p_value))
}) %>% bind_rows()

# 88% of the time, the individuals have significantly higher genetic diversity at adaptive SNPs than neutral SNPs
pvals %>%
  mutate_at(vars(everything()), ~p.adjust(.x, method = "fdr")) %>%
  mutate_at(vars(everything()), ~(.x < 0.05)) %>%
  summarise_all(mean)

# Calculate correlation between adaptive and neutral heterozygosity
t_het_neutral <- list_transpose(het_neutral2)
rvals <- map_dbl(1:length(t_het_neutral), ~{
  cor(het_adaptive, t_het_neutral[[.x]], use = "pairwise.complete.obs")
}) 
# The mean correlation between adaptive and neutral heterozygosity is 0.67
mean(rvals)

mean_neutral <- map_dbl(het_neutral, ~mean(.x, na.rm = TRUE))
mean_neutral2 <- map_dbl(het_neutral2, ~mean(.x, na.rm = TRUE))
gg_df <-
  coords %>%
  mutate(neutral1 = mean_neutral, neutral2 = mean_neutral2, adaptive = het_adaptive) %>%
  pivot_longer(cols = c(neutral1, neutral2), names_to = "maf", values_to = "neutral") %>%
  mutate(maf = ifelse(maf == "neutral1", "no filter", "maf filter"))

pdf(here("analysis", "gea", "het_adaptive_vs_neutral2.pdf"), width = 5, height = 5)
ggplot(gg_df, aes(x = neutral, y = adaptive)) +
  geom_point(aes(color = maf)) +
  geom_abline(intercept = 0, slope = 1, color = "black", lty = "dashed") +
  geom_smooth(method = "lm", se = TRUE, aes(color = maf)) +
  labs(color = "Filter") +
  labs(x = "Neutral SNP Heterozygosity", y = "Adaptive SNP Heterozygosity") +
  theme_classic()

gg_long <- pivot_longer(gg_df, cols = c(neutral, adaptive), names_to = "type", values_to = "het")
ggplot(gg_long) +
  geom_sf(aes(color = het)) +
  scale_color_viridis_c(option = "turbo") +
  theme_void() +
  facet_wrap(~type + maf) +
  ggtitle("Heterozygosity at Adaptive and Neutral SNPs (Raw)")

gg_scaled <-
  gg_df %>%
  mutate(neutral = scale(neutral), adaptive = scale(adaptive)) %>%
  pivot_longer(cols = c(neutral, adaptive), names_to = "type", values_to = "het")
ggplot(gg_scaled) +
  geom_sf(aes(color = het)) +
  scale_color_viridis_c(option = "turbo") +
  theme_void() +
  facet_wrap(~type + maf) +
  ggtitle("Heterozygosity at Adaptive and Neutral SNPs (Scaled)")

# Plot across space
dev.off()

# TEST 2: Make frequencies of SNPs in distribution the same as those of adaptive SNPs
freq_adaptive <- apply(dos[,rda_snps], 2, function(x) mean(x, na.rm = TRUE)/2)
dos_neutral <- dos[,which(!(colnames(dos) %in% rda_snps))]
freq_neutral <- apply(dos_neutral, 2, function(x) mean(x, na.rm = TRUE)/2)

# PROBLEM WITH THIS: might just be selecting the same SNPs multiple times
picks <- map(rda_snps, \(snp){
    freq  <- round(freq_adaptive[snp], 1)
    neutral_options <- freq_neutral[freq_neutral == freq]
    if (length(neutral_options) == 0) return(NA)
    return(names(neutral_options))
  })

set.seed(354)
het_neutral3 <- map(1:100, ~{
  idx <- map_chr(picks, ~sample(.x, 1))
  warning("No neutral options for for ", round(mean(is.na(idx)), 2)*100, "% of SNPs")
  dos_to_het(dos_neutral[, na.omit(idx)])
}, .progress = TRUE) %>% 
list_transpose()

set.seed(354)
freq_check <- map(1:100, ~{
  idx <- map_chr(picks, ~sample(.x, 1))
  warning("No neutral options for for ", round(mean(is.na(idx)), 2)*100, "% of SNPs")
  apply(dos_neutral[, na.omit(idx)], 2, function(x) mean(x, na.rm = TRUE)/2)
}, .progress = TRUE) 
all(freq_check[[1]] == freq_check[[2]])
!all(names(freq_check[[1]]) == names(freq_check[[2]]))
!all(dos[1, names(freq_check[[1]])] == dos[1, names(freq_check[[2]])])

# For each individual, calculate a p-value
pvals <- map(names(het_adaptive), ~{
  value <- het_adaptive[.x]
  p_greater <- mean(het_neutral3[[.x]] >= value)
  p_less <- mean(het_neutral3[[.x]] <= value)
  p_value <- 2 * min(p_greater, p_less)
  return(c(greater = p_greater, less = p_less, both = p_value))
}) %>% bind_rows()

# 88% of the time, the individuals have significantly higher genetic diversity at adaptive SNPs than neutral SNPs
pvals %>%
  mutate_at(vars(everything()), ~p.adjust(.x, method = "fdr")) %>%
  mutate_at(vars(everything()), ~(.x < 0.05)) %>%
  summarise_all(mean)

# Calculate correlation between adaptive and neutral heterozygosity
t_het_neutral <- list_transpose(het_neutral3)
rvals <- map_dbl(1:length(t_het_neutral), ~{
  cor(het_adaptive, t_het_neutral[[.x]], use = "pairwise.complete.obs")
}) 
# The mean correlation between adaptive and neutral heterozygosity is 0.67
mean(rvals)

mean_neutral <- map_dbl(het_neutral, ~mean(.x, na.rm = TRUE))
mean_neutral2 <- map_dbl(het_neutral2, ~mean(.x, na.rm = TRUE))
mean_neutral3 <- map_dbl(het_neutral3, ~mean(.x, na.rm = TRUE))


gg_df <-
  coords %>%
  mutate(neutral1 = mean_neutral, neutral2 = mean_neutral2, neutral3 = mean_neutral3, adaptive = het_adaptive) %>%
  pivot_longer(cols = c(neutral1, neutral2, neutral3), names_to = "maf", values_to = "neutral") %>%
  mutate(maf = case_when(
    maf == "neutral1" ~ "no filter",
    maf == "neutral2" ~ "maf filter",
    maf == "neutral3" ~ "freq distribution"
    )
  )

pdf(here("analysis", "gea", "het_adaptive_vs_neutral3.pdf"), width = 6.5, height = 5)
ggplot(gg_df, aes(x = neutral, y = adaptive)) +
  geom_point(aes(color = maf), alpha = 0.5, cex = 0.7, cex = 16) +
  geom_abline(intercept = 0, slope = 1, color = "black", lty = "dashed") +
  geom_smooth(method = "lm", se = TRUE, aes(color = maf)) +
  labs(color = "Filter") +
  labs(x = "Neutral SNP Heterozygosity", y = "Adaptive SNP Heterozygosity") +
  theme_classic()

gg_long <- pivot_longer(gg_df, cols = c(neutral, adaptive), names_to = "type", values_to = "het")

ggplot(gg_long) +
  geom_density(aes(x = het, fill = maf, color = type), alpha = 0.5) +
  theme_classic() +
  ggtitle("Heterozygosity at Adaptive and Neutral SNPs")

ggplot(gg_long) +
  geom_sf(aes(color = het)) +
  scale_color_viridis_c(option = "turbo") +
  theme_void() +
  facet_wrap(~type + maf) +
  ggtitle("Heterozygosity at Adaptive and Neutral SNPs (Raw)")

gg_scaled <-
  gg_df %>%
  mutate(neutral = scale(neutral), adaptive = scale(adaptive)) %>%
  pivot_longer(cols = c(neutral, adaptive), names_to = "type", values_to = "het")
ggplot(gg_scaled) +
  geom_sf(aes(color = het)) +
  scale_color_viridis_c(option = "turbo") +
  theme_void() +
  facet_wrap(~type + maf) +
  ggtitle("Heterozygosity at Adaptive and Neutral SNPs (Scaled)")

dev.off()



# IDEA:
1. Identify SNPs associated with adaptation to warmer temperatures
2. Predict where allele frequencies would need to increase to adapt to warmer temperatures
3. Compare the genetic diversity at these SNPs to the genetic diversity at neutral SNPs

# Compare with whole genome heterozygosity
source(here("analysis", "genetic_diversity", "genetic_diversity.R"))
het <- get_het()
gg_df <-
  coords %>% 
  left_join(het) %>%
  mutate(neutral = scale(Ho), adaptive = scale(het_adaptive))

pdf(here("analysis", "gea", "het_adaptive_vs_neutral4.pdf"), width = 6, height = 5)
ggplot(gg_df, aes(x = neutral, y = adaptive)) +
  geom_point(alpha = 0.5) +
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