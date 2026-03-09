library(tidyverse)
library(here)
source(here("general_functions.R"))
pop_df <- get_pops()
# popmap.txt
Ind1 Ind1 pop1
Ind2 Ind2 pop1
Ind3 Ind3 pop2
Ind4 Ind4 pop2

pop_plink <- 
  pop_df %>% 
  select(IID = SampleID, POP = cluster) %>% 
  mutate(FID = 0) %>%
  select(FID, IID, POP)

head(pop_plink)

write_delim(pop_plink, 
            delim = " ", 
            file = here("analysis", "check_nonsyn", "outputs", "popmap.txt"), 
            col_names = FALSE)
