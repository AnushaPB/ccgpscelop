# This is a test for myself to make sure that the Southern CA populations are not driving the patterns we see in the genetic diversity analyses. I will run the same models as before but excluding the Southern CA populations (clusters 6 and 8) and see if the results are consistent.
library(tidyverse)
library(here)
source(here("general_functions.R"))
model_df <- read_csv(here("analysis", "genetic_diversity", "outputs", "model_df.csv"))
head(model_df)
pop_df <- get_pops()
model_df <- left_join(model_df, pop_df)

no_socal <- model_df %>% filter(!cluster %in% c(6, 8))
lm(Ho ~ csi_past + bio1 + tmean_dif + Q + fire_frq + fire_recent, data = no_socal) %>% summary()
unique(no_socal$cluster)

cluster_cols <- viridis::viridis(n = length(unique(model_df$cluster)), option = "turbo")
names(cluster_cols) <- sort(unique(model_df$cluster))
quick_plot <- function(x) {
  plt1 <- 
    ggplot(data = model_df, aes_string(x = x, y = "Ho")) +
    geom_point(aes(color = factor(cluster))) +
    geom_smooth(method = "lm", se = FALSE) +
    theme_classic() +
    labs(color = "Cluster") +
    scale_color_manual(values = cluster_cols) +
    ggpubr::stat_cor()
  
  plt2 <-
    ggplot(data = no_socal, aes_string(x = x, y = "Ho")) +
    geom_point(aes(color = factor(cluster))) +
    geom_smooth(method = "lm", se = FALSE) +
    theme_classic() +
    labs(color = "Cluster") +
    scale_color_manual(values = cluster_cols) +
    ggpubr::stat_cor()

  cowplot::plot_grid(plt1, plt2, labels = c("A. All populations", "B. Excluding Southern CA"), ncol = 2)
}

library(cowplot)
png(here("analysis", "genetic_diversity", "plots", "southern_ca_test.png"), width = 10, height = 20, res = 300, units = "in")
plot_grid(quick_plot("csi_past"), quick_plot("bio1"), quick_plot("tmean_dif"), quick_plot("fire_frq"), quick_plot("fire_recent"), ncol = 1)
dev.off()

mod <- lm(Ho ~ csi_past + bio1 + tmean_dif + evt + fire_frq * fire_recent + cluster, data = model_df) %>% summary()
broom::tidy(mod) %>% filter(p.value < 0.05, !term %in% c("(Intercept)"))


mod <- lm(Ho ~ csi_past + bio1 + tmean_dif + evt + fire_frq * fire_recent + cluster * tmean_dif, data = model_df) %>% summary()
broom::tidy(mod) %>% filter(p.value < 0.05, !term %in% c("(Intercept)")) 

broom::tidy(mod) %>% filter(p.value < 0.05, !term %in% c("(Intercept)"))  %>% filter(str_detect(term, "tmean_dif"))

# plot tmean_dif vs Ho by cluster
pdf(here("tmean_dif_by_cluster.pdf"), width = 10, height = 10)
plt1 <-
  ggplot(data = model_df, aes(x = tmean_dif, y = Ho)) +
  geom_point(aes(color = factor(cluster))) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(color = "Cluster") +
  scale_color_manual(values = cluster_cols) +
  ggpubr::stat_cor() +
  facet_wrap(~ cluster, scales = "free")

# Map of clusters
plt2 <- 
  ggplot(data = model_df, aes(x = x, y = y)) +
  geom_point(aes(color = factor(cluster))) +
  theme_classic() +
  labs(color = "Cluster") +
  scale_color_manual(values = cluster_cols)

plt2
plt1
dev.off()
