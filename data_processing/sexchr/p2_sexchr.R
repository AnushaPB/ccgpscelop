library(here)
library(tidyverse)
source(here("general_functions.R"))

coords <- get_coords()
original_meta <- 
  read_csv(here("data", "sceloporus_metadata.csv")) %>%
  select(SampleID, MinicoreID, sex)

museum <- 
  read_csv(here("data_processing", "sexchr", "all_museum.csv")) %>% 
  select(SampleID, museum_sex = sex) %>% 
  filter(SampleID %in% original_meta$SampleID)
table(museum$museum_sex)

meta <- 
  original_meta %>% 
  mutate(sex = ifelse(sex == "NaN", NA, sex)) %>%
  left_join(museum) %>%
  mutate(museum_sex = ifelse(museum_sex == "undetermined", NA, museum_sex)) %>%
  mutate(sex = ifelse(is.na(sex), museum_sex, sex))
table(meta$sex)

males <- meta %>% filter(sex == "male") %>% filter(MinicoreID != "NaN") %>% pull(MinicoreID)
females <- meta %>% filter(sex == "female") %>% filter(MinicoreID != "NaN") %>% pull(MinicoreID)
writeLines(males, here("data_processing", "sexchr", "males.txt"))
writeLines(females, here("data_processing", "sexchr", "females.txt"))

# Depth
sex1_depth <- read_table(here("data_processing", "sexchr", "sex_linked_1_depth.idepth")) %>% mutate(chr = "sex_linked_1")
sex2_depth <- read_table(here("data_processing", "sexchr", "sex_linked_2_depth.idepth")) %>% mutate(chr = "sex_linked_2")
sex_depth <- bind_rows(sex1_depth, sex2_depth) %>% rename(MinicoreID = INDV) %>% select(-N_SITES)

# Missingness 
sex1_miss <- read_table(here("data_processing", "sexchr", "sex_linked_1.imiss")) %>% mutate(chr = "sex_linked_1")
sex2_miss <- read_table(here("data_processing", "sexchr", "sex_linked_2.imiss")) %>% mutate(chr = "sex_linked_2")
sex_miss <- bind_rows(sex1_miss, sex2_miss) %>% rename(MinicoreID = INDV) %>% select(MinicoreID, chr, F_MISS)

# Heterozygosity
sex1_het <- read_table(here("data_processing", "sexchr", "sex_linked_1.het")) %>% mutate(chr = "sex_linked_1")
sex2_het <- read_table(here("data_processing", "sexchr", "sex_linked_2.het")) %>% mutate(chr = "sex_linked_2")
sex_het <- 
  bind_rows(sex1_het, sex2_het) %>% 
  rename(MinicoreID = IID) %>% 
  dplyr::select(-FID) %>% 
  mutate(Ho = (`N(NM)` - `O(HOM)`) / `N(NM)`) 

sex_df <- left_join(sex_het, sex_depth) %>% left_join(sex_miss) %>% left_join(meta)

table(sex_df$sex)



plt1 <- 
  ggplot(sex_df) +
  geom_point(aes(x = Ho, y = MEAN_DEPTH), alpha = 0.8, col = "lightgray") +
  geom_point(data = drop_na(sex_df, sex), aes(x = Ho, y = MEAN_DEPTH, col = sex)) +
  theme_minimal() +
  facet_wrap(~chr, nrow=1) +
  labs(x = "Heterozygosity", y = "Mean Depth") 

plt2 <-
  ggplot(sex_df) +
  geom_point(aes(x = Ho, y = F_MISS), alpha = 0.8, col = "lightgray") +
  geom_point(data = drop_na(sex_df, sex), aes(x = Ho, y = F_MISS, col = sex)) +
  theme_minimal() +
  facet_wrap(~chr, nrow=1) +
  labs(x = "Heterozygosity", y = "Missingness") 

plt3 <- 
  ggplot(sex_df) +
  geom_point(aes(x = MEAN_DEPTH, y = F_MISS), alpha = 0.8, col = "lightgray") +
  geom_point(data = drop_na(sex_df, sex), aes(x = MEAN_DEPTH, y = F_MISS, col = sex)) +
  theme_minimal() +
  facet_wrap(~chr, nrow=1) +
  labs(x = "Mean Depth", y = "Missingness") 

pdf(here("data_processing", "sexchr", "sexchr_plot.pdf"), width=7, height=9)
cowplot::plot_grid(plt3, plt2, plt1, nrow=3)
dev.off()

sex_df %>% drop_na(sex) %>% group_by(chr, sex) %>% summarize_at(c("Ho", "MEAN_DEPTH", "F_MISS"), mean)

# Summary:
# sex_linked_1 - heterozygosity higher in females than males, depth higher in females than males, and missingness lower in females than males -> probably X chromosome
# sex_linked_2 - heterozygosity higher in males than females, depth higher in males than females, and missingness lower in males than females -> probably Y chromosome
