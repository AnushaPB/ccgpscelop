library(here)
library(tidyverse)
library(sf)
devtools::load_all()

# check depth and missingness, averaged across all SNPS

# PART 1: CHECKING FULL SNP DATASET (POST DEPTH AND MAF FILTER) -------------------------------------------------
# get data
path <- here("58-Sceloporus")
depth <- read_table(here(path, "sample_depth_info.idepth"))
miss <- read_table(here(path, "sample_missing_info.imiss"))

coords <- read_table(here("data", "58-Sceloporus.coords.txt"), col_names = FALSE)
colnames(coords) <- c("INDV", "x", "y")
coords <- sf::st_as_sf(coords, coords = c("x", "y"))

df <-
  left_join(depth, miss) %>%
  right_join(coords)

# plotting depth/missingness
ggplot(df) +
  geom_point(aes(x = F_MISS, y = MEAN_DEPTH)) +
  #geom_vline(xintercept = 0.6, lty = "dashed") +
  theme_classic()

# removing all individuals with greater than 20% missing data
df %>% filter(F_MISS > 0.20) %>% pull(INDV)

# plot where the individuals are
ggplot(df) +
  geom_sf(aes(geometry = geometry)) +
  geom_sf(data = filter(df, F_MISS > 0.20), aes(geometry = geometry), col = "red")

# repeat for vira:
ggplot(vira_df) +
  geom_point(aes(x = F_MISS, y = MEAN_DEPTH, col = group)) +
  geom_vline(xintercept = 0.2, lty = "dashed") +
  theme_classic()

# one individual with missingness greater than 0.20:"Ralllimi_CCGPMC042_HSU9048_RL"
vira_miss %>% filter(F_MISS > 0.20) %>% pull(INDV)

ggplot(vira_df) +
  geom_sf(aes(geometry = geometry)) +
  geom_sf(data = filter(vira_df, F_MISS > 0.20), aes(geometry = geometry), col = "red")

# PART 2: CHECK FILTERED DATA (maf > 0.05, DP > 5, DP < 50)
path <- here("data", "ccgp_data", "QC")
blra_miss <- read.table(here(path, "Laterallus_postfilter_sample_missing_info.imiss"), header = 1)
vira_miss <- read.table(here(path, "Rallus_postfilter_sample_missing_info.imiss"), header = 1)

coords <- get_coords()

blra_df <-
  blra_miss %>%
  mutate(SampleID = recode_blra(INDV)) %>%
  right_join(coords$blra)

vira_df <-
  vira_miss %>%
  mutate(SampleID = recode_vira(INDV)) %>%
  right_join(coords$vira)

# plotting depth/missingness
ggplot(blra_df) +
  geom_histogram(aes(x = F_MISS, fill = group)) +
  geom_vline(xintercept = 0.2, lty = "dashed") +
  theme_classic()

# multiple individuals with missingness greater than 0.60:
#[1] "Latejama_CCGPMC040_92168903_LJ" "Latejama_CCGPMC040_92168908_LJ"
#[3] "Latejama_CCGPMC040_92168913_LJ" "Latejama_CCGPMC040_92168914_LJ"
#[5] "Latejama_CCGPMC040_92168915_LJ" "Latejama_CCGPMC040_92168921_LJ"
#[7] "Latejama_CCGPMC040_92168923_LJ" "LatJam_CCGPMC038_92168988_LJ"  
# removing all individuals with greater than 60% missing data
blra_miss %>% filter(F_MISS > 0.60) %>% pull(INDV)

ggplot(blra_df) +
  geom_sf(aes(geometry = geometry)) +
  geom_sf(data = filter(blra_df, F_MISS > 0.60), aes(geometry = geometry), col = "red")

# repeat for vira
ggplot(vira_df) +
  geom_histogram(aes(x = F_MISS, fill = group)) +
  geom_vline(xintercept = 0.2, lty = "dashed") +
  theme_classic()

# two individual with missingness greater than 0.60: "RalLim_CCGPMC031_171304172_RL" "RalLim_CCGPMC031_171304195_RL"
vira_miss %>% filter(F_MISS > 0.60) %>% pull(INDV)

ggplot(vira_df) +
  geom_sf(aes(geometry = geometry)) +
  geom_sf(data = filter(vira_df, F_MISS > 0.60), aes(geometry = geometry), col = "red")

# PART 3: REMOVE INDIVIDUALS ------------------------------------------------------------
# remove individuals and write out a text file for each population
# matched by these steps in data_processing.sh:
#bcftools view -s ^RalLim_171304107_RL,RalLim_171304111_RL,Ralllimi_CCGPMC042_HSU9048_RL 9-Rallus_snpsonly.recode.vcf -Oz -o 9-Rallus_snpsonly_rmsamp.vcf
#bcftools view -s ^LatJam_Z18-430_LJ,LatJam_Z20-714_LJ,Latejama_CCGPMC041_242161856_LJ 9-Laterallus_snpsonly.recode.vcf -Oz -o 9-Laterallus_snpsonly_rmsamp.vcf
vira_rmsamp <-
  vira_df %>%
  # from PART 1: removing individuals with > 20% missing data in pre-filter data (SNPs only)
  filter(!(INDV %in% c("RalLim_171304107_RL", "RalLim_171304111_RL", "CCGPMC042_HSU9048_RL"))) %>%
  # from PART 2: removing individuals with > 60% missing data in post-filter data
  filter(!(INDV %in% c("RalLim_CCGPMC031_171304172_RL", "RalLim_CCGPMC031_171304195_RL")))
  
vira_sf <- 
  vira_rmsamp %>%
  filter(group == "Sierra Foothills") %>%
  pull(INDV)

vira_ba <-
  vira_rmsamp %>%
  filter(group == "Bay Area") %>%
  pull(INDV)

blra_rmsamp <-
  blra_df %>%
  # from PART 1: removing individuals with > 20% missing data in pre-filter data (SNPs only)
  filter(!(INDV %in% c("LatJam_Z18-430_LJ", "LatJam_Z20-714_LJ", "Latejama_CCGPMC041_242161856_LJ"))) %>%
  # from PART 2: removing individuals with > 60% missing data in post-filter data
  filter(!(INDV %in% c(
    "Latejama_CCGPMC040_92168903_LJ", 
    "Latejama_CCGPMC040_92168908_LJ", 
    "Latejama_CCGPMC040_92168913_LJ", 
    "Latejama_CCGPMC040_92168914_LJ",
    "Latejama_CCGPMC040_92168915_LJ",
    "Latejama_CCGPMC040_92168921_LJ",
    "Latejama_CCGPMC040_92168923_LJ",
    "LatJam_CCGPMC038_92168988_LJ"
    ))) 

blra_sf <-
  blra_rmsamp %>%
  filter(group == "Sierra Foothills") %>%
  pull(INDV)

blra_ba <-
  blra_rmsamp %>%
  filter(group == "Bay Area") %>%
  pull(INDV)

blra_sc <-
  blra_rmsamp %>%
  filter(group == "S. CA") %>%
  pull(INDV)

message("writing out _samples.txt files")
writeLines(vira_sf, here("data", "ccgp_data", "QC", "vira_sf_samples.txt"))
writeLines(vira_ba, here("data", "ccgp_data", "QC", "vira_ba_samples.txt"))
writeLines(blra_sf, here("data", "ccgp_data", "QC", "blra_sf_samples.txt"))
writeLines(blra_ba, here("data", "ccgp_data", "QC", "blra_ba_samples.txt"))
writeLines(blra_sc, here("data", "ccgp_data", "QC", "blra_sc_samples.txt"))
