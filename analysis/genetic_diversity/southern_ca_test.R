library(tidyverse)
source(here("general_functions.R"))
model_df <- read_csv(here("analysis", "genetic_diversity", "outputs", "model_df.csv"))
head(model_df)
pop_df <- get_pops()
model_df <- left_join(model_df, pop_df)

no_socal <- model_df %>% filter(!cluster %in% c(6, 8))
lm(Ho ~ csi_past + bio1 + tmean_dif + Q + log_fri + pfs_surfac, data = no_socal) %>% summary()
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

png(here("analysis", "genetic_diversity", "plots", "southern_ca_test.png"), width = 10, height = 20, res = 300, units = "in")
plot_grid(quick_plot("csi_past"), quick_plot("bio1"), quick_plot("tmean_dif"), quick_plot("log_fri"), quick_plot("pfs_surfac"), ncol = 1)
dev.off()
