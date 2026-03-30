library(tidyverse)
library(sf)
library(here)
source(here("general_functions.R"))

get_coords() %>%
  select(SampleID) %>%
  write_tsv(here("data", "final_sampleids.txt"), col_names = FALSE)
