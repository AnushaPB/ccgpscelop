library(here)
library(tidyverse)
library(vcfR)
library(wingen)
library(terra)
library(viridis)
library(tmap)
# DISTINCT SIM 2
filepath <- here("gnx", "GNX_mod-distinct_gradient_sim", "it-0", "spp-spp_0")
gsd <- get_gsd(here(filepath, "mod-distinct_gradient_sim_it-0_t-1000_spp-spp_0.csv"))
vcf <- read.vcfR(here(filepath, "mod-distinct_gradient_sim_it-0_t-1000_spp-spp_0.vcf"))
dos <- vcf_to_dosage(vcf)
freqs <- dos/2

df <- 
bind_cols(gsd, dos[,1:8]) %>%
pivot_longer(starts_with("0_"))

ggplot(gsd) +
  geom_point(aes(x = x, y = y, col = z1)) +
  coord_fixed()

ggplot(gsd) +
  geom_point(aes(x = x, y = y, col = envlyr, pch = factor(klyr))) +
  coord_fixed()

ggplot(df) +
  geom_point(aes(x = x, y = y, col = value)) +
  facet_wrap(~name, nrow = 2) +
  coord_fixed()

het <- freqs[s,]
het[het == 1] <- 0
het[het == 0] <- 0
het[het == 0.5] <- 1
het <- apply(het, 1, mean)

subgsd <- gsd[s,]
subgsd$het <- het
subgsd$adaptive <- distinct_ind(freqs[s, 1:4], rare_allele = TRUE)[[1]]
subgsd$neutral <- distinct_ind(freqs[s, 5:8], rare_allele = TRUE)[[1]]

ggdf <- subgsd %>% pivot_longer(c("adaptive", "neutral", "resid"))
ggplot(ggdf) + geom_point(aes(x = x, y = y, col = value)) + scale_color_viridis_c(option = "turbo") + facet_wrap(~name)

ggplot(subgsd) + geom_point(aes(x = x, y = y, col = het)) + scale_color_viridis_c(option = "turbo")

mod <- lm(adaptive ~ neutral, subgsd)
subgsd$resid1 <- resid(mod)
mod <- lm(adaptive ~ het + neutral, subgsd)
subgsd$resid <- resid(mod)
