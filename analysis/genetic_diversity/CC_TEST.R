library(tidyverse)

df <- read.csv("analysis/genetic_diversity/outputs/model_df.csv")
pops <- get_pops()
df <- df %>% left_join(pops) %>% st_as_sf(coords = c("x", "y"), crs = 4326) %>% st_transform(3310)

mod <- lm(Ho ~ tmean_dif, data = df)
summary(mod)

pdf(here("analysis/genetic_diversity/plots/Ho_tmean_dif_by_pop.pdf"), height = 6, width = 8)

ggplot(df, aes(x = tmean_dif, y = Ho, color = cluster)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Ho vs. Tmean Dif by Population", x = "Tmean Dif", y = "Ho") +
  theme_classic()

dev.off()

library(ggplot2)
library(dplyr)
library(broom)
library(purrr)
library(ggpubr)
library(here)

# Fit linear models by cluster and extract significant ones
models <- df %>%
  group_by(cluster) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(Ho ~ tmean_dif, data = .x)),
    tidied = map(model, tidy)
  ) %>%
  unnest(tidied) %>%
  filter(term == "tmean_dif", p.value < 0.1) %>%
  pull(cluster)

# Plot
library(ggplot2)
library(dplyr)
library(emmeans)
library(cowplot)
library(here)

# Step 1: Fit model with interaction
mod <- lm(Ho ~ tmean_dif * cluster, data = df)

# Step 2: Post hoc test - which clusters have a significant slope for tmean_dif?
library(emmeans)
em_trends <- emtrends(mod, var = "tmean_dif", specs = "cluster")
sig_clusters <- 
  summary(em_trends, infer = c(TRUE, TRUE)) %>%
  filter(p.value < 0.1) %>%
  pull(cluster)

# Step 3: Plots
plt1 <- ggplot(df, aes(x = tmean_dif, y = Ho, color = cluster)) +
  geom_point() +
  geom_smooth(
    method = "lm", se = FALSE,
    data = df %>% filter(cluster %in% sig_clusters)
  ) +
  labs(title = "Ho vs. Tmean Dif by Population (Significant Slopes Only)",
       x = "Tmean Dif", y = "Ho") +
  theme_classic() +
  theme(legend.position = "right")

plt2 <- ggplot(df) + 
  geom_sf(aes(color = cluster)) +
  theme_void()

# Step 4: Save side-by-side plots
pdf(here("analysis/genetic_diversity/plots/Ho_tmean_dif_by_pop.pdf"), height = 6, width = 10)
cowplot::plot_grid(plt1, plt2, nrow = 1)
dev.off()

mod <- lm(Ho ~ tmean_dif * bio1 + csi_past + Q + NDVI + gHM, data = df)

mod2 <- lm(Ho ~ tmean_dif + bio1 + csi_past + Q + NDVI + gHM, data = df)

summary(mod)
summary(mod2)

AIC(mod, mod2)

# Create a grid of predictor values
grid <- expand.grid(
  tmean_dif = seq(min(df$tmean_dif), max(df$tmean_dif), length.out = 100),
  bio1 = quantile(df$bio1, probs = seq(0, 1, 0.1)), 
  csi_past = mean(df$csi_past)
)

# Predict Ho from model
grid$Ho <- predict(mod, newdata = grid)

# Plot
pdf(here("analysis/genetic_diversity/plots/Ho_tmean_dif_bio1.pdf"), height = 4, width = 5.5)

ggplot(grid, aes(x = tmean_dif, y = Ho, color = bio1, group = factor(bio1))) +
  geom_line(size = 1) +
  labs(
    x = "Climate Change",
    y = "Predicted Heterozygosity (Ho)",
    color = "Contemporary\ntemperature"
  ) +
  scale_color_gradientn(colors = rev(MetBrewer::met.brewer("Hiroshige", type = "continuous"))) +
  theme_classic()
dev.off()

library(dplyr)
library(ggplot2)
library(here)library(dplyr)
library(tibble)
library(ggplot2)
library(here)

# 1. Compute mean bio1 per cluster
cluster_bio1 <- df %>%
  group_by(cluster) %>%
  summarize(mean_bio1 = mean(bio1, na.rm = TRUE), .groups = "drop")

# 2. Create prediction grid for each cluster
grid <- cluster_bio1 %>%
  rowwise() %>%
  mutate(pred_data = list(
    tibble(
      tmean_dif = seq(min(df$tmean_dif), max(df$tmean_dif), length.out = 100),
      bio1 = mean_bio1,
      csi_past = mean(df$csi_past, na.rm = TRUE),
      cluster = cluster
    )
  )) %>%
  pull(pred_data) %>%
  bind_rows()

# 3. Predict Ho from model
grid$Ho <- predict(mod, newdata = grid)

# 4. Plot
pdf(here("analysis/genetic_diversity/plots/Ho_tmean_dif_bio1_by_cluster.pdf"), height = 6, width = 8)

ggplot(grid, aes(x = tmean_dif, y = Ho, color = cluster)) +
  geom_line(size = 1) +
  labs(
    title = "Predicted Heterozygosity vs. Climate Change",
    subtitle = "Each line reflects a cluster's average temperature (bio1)",
    x = "Climate Change (tmean_dif)",
    y = "Predicted Heterozygosity (Ho)",
    color = "Cluster"
  ) +
  theme_classic()

dev.off()