library(tidyverse)
library(here)

chr_key <- 
  read_csv(here("data_processing", "chromosemble", "outputs", "chromosome_labels_p.csv")) %>%
  mutate(
    chr_names_v1.0 = case_when(scaffold == "SCAF_7" ~ "sex_linked_1", scaffold == "SCAF_8" ~ "sex_linked_2", TRUE ~ chr),
    chr_names_v2.0 = case_when(scaffold == "SCAF_7" ~ "X", scaffold == "SCAF_8" ~ "Y", TRUE ~ chr)
  ) %>%
  select(-chr)

write_csv(chr_key, here("data_processing", "sexchr", "outputs", "sex_chromosome_key.csv"))
  


