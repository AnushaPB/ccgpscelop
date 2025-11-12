library(tidyverse)
library(here)
source(here("analysis", "genetic_diversity", "functions_genetic_diversity.R"))

gea <- get_het("nonsyn.het") %>% select(SampleID, gea_Ho = Ho)
mod_df <- read_csv(here("analysis", "genetic_diversity", "outputs", "model_df.csv")) 
df <- left_join(mod_df, gea)

pdf(here("TEMP.pdf"))
ggplot(df, aes(x = Ho, y = gea_Ho)) +
  geom_point(aes(col = bio1)) +
  geom_smooth(method = "lm") +
  theme_classic() +
  scale_color_viridis_c(option = "turbo") +
  ggpubr::stat_cor()
dev.off()

summary(lm(gea_Ho ~ Ho + bio1 + tmean_dif + csi_past + NDVI + glacier + gHM + Q, data = df))
