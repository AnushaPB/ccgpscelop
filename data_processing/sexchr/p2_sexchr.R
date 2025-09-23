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

# Depth across entire genome
genomewide <- 
  read_table(here("data_processing", "sexchr", "outputs", "genomewide.idepth")) %>% 
  rename(MinicoreID = INDV, genomewide_depth = MEAN_DEPTH) %>% 
  select(-N_SITES)

# Samples in the same order as the DP columns
samples <- scan(here("data_processing", "sexchr", "outputs", "samples.txt"), what = "")
col_names <- c("CHROM", "POS", samples)

sex_depth_ind <- read_tsv(
  here("data_processing", "sexchr", "outputs", "sex_linked.dp.tsv"),
  col_names = col_names,
  na = c(".", "NA", ""),                # <- treat "." as NA
  col_types = cols(
    CHROM = col_character(),
    POS   = col_integer(),
    .default = col_double()             # all DP columns numeric
  )
) %>%
pivot_longer(
  cols = all_of(samples),
  names_to = "SampleID",
  values_to = "DP"
) 

sex_depth_ind_window <- 
  sex_depth_ind %>%
  mutate(window = (POS %/% 100000) * 100000) %>%
  group_by(SampleID, window, CHROM) %>%
  summarize(mean_DP = mean(DP, na.rm=TRUE)) %>%
  ungroup() 
  
window_depth_df <-
  sex_depth_ind_window %>%
  rename(MinicoreID = SampleID) %>%
  left_join(meta, by = "MinicoreID") %>%
  left_join(genomewide, by = c("MinicoreID")) %>%
  mutate(scaled_depth = mean_DP / genomewide_depth)

averaged_by_sex <- 
  window_depth_df %>%
  drop_na(sex) %>%
  group_by(window, CHROM, sex) %>%
  summarize(scaled_depth = mean(scaled_depth, na.rm=TRUE)) 

pdf(here("data_processing", "sexchr", "sexchr_window_depth.pdf"), width=8, height=4)
ggplot(window_depth_df) +
  geom_line(
    data = filter(window_depth_df, sex == "male"), # "female" before "male"
    aes(x = window, y = scaled_depth, col = sex, group = MinicoreID),
    alpha = 0.05) + 
  geom_line(
    data = filter(window_depth_df, sex == "female"), # "female" before "male"
    aes(x = window, y = scaled_depth, col = sex, group = MinicoreID),
    alpha = 0.05
  ) +
  geom_line(data = averaged_by_sex, aes(x = window, y = scaled_depth, col = sex), linewidth = 1) +
  theme_minimal() +
  labs(x = "Genomic Window", y = "Scaled Mean Depth (100kb Windows)") +
  facet_wrap(~CHROM, ncol = 1)
dev.off()
