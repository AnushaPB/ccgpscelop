library(here)
library(tidyverse)
library(sf)
devtools::load_all()
source(here("general_functions.R"))

path <- here("data", "ccgp_data")
coords <- get_coords() %>% rename(INDV = SampleID)

# check depth and missingness, averaged across all SNPS
ids <- read_table(here(path, "58-Sceloporus_annotated_pruned_0.6.fam"), col_names = FALSE)$X2
depth <- read_table(here(path, "58-Sceloporus.idepth"))
miss <- read_table(here(path, "58-Sceloporus.imiss"))


df <-
  left_join(depth, miss) %>%
  right_join(coords) %>%
  # Remove any individuals already dropped from clean snps
  filter(INDV %in% ids) 

bad <-
  df %>%
  filter(F_MISS > 0.2 | MEAN_DEPTH < 8)

pdf(here("data_processing", "depth_missingness.pdf"))
# plotting depth/missingness
ggplot(df, aes(x = F_MISS, y = MEAN_DEPTH)) +
  geom_point() +
  geom_point(data = bad, color = "red") +
  ggrepel::geom_text_repel(data = bad, aes(label = INDV), color = "red") +
  theme_classic()
dev.off()

# JEM16−001−S14_na is the one major problem individual (high missingness, low depth)
# Scelocci_MVZ_265296 and Scelocci_CAS219619 also have low depth
bad %>% select(INDV, F_MISS, MEAN_DEPTH)
