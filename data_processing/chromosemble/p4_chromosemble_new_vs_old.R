# Satsuma alignment -------------------------------------------------------
# Modified from code provided Erik Enbody 
library(tidyverse)
library(here)

# Load chain file
chain.ens <- read_table(here("data_processing", "chromosemble2","new_vs_old", "xcorr_first", "satsuma_summary.chained.out"), col_names = FALSE)
names(chain.ens) <- c("qName", "qStart", "qStop", "refName", "refStart", "refStop", "V7", "dir")
head(chain.ens)

# Format data
chain.ens.out <- 
  chain.ens %>% 
  # Remove empty alignments
  filter(refName!="") %>% 
  # Calculate length of alignment
  mutate(qLength = qStop - qStart)

# Get orientation by choosing the most common orientation for each scaffold
df.orientation <- 
  chain.ens.out %>% 
  group_by(qName, dir) %>% 
  # Count how many alignments in each direction
  summarise(n.dir = n()) %>% 
  group_by(qName) %>% 
  # Selection most common direction
  filter(n.dir == max(n.dir)) %>% 
  group_by(qName) %>% 
  # For those with few alignments, sometimes two orientations have same number of matches. This just selects 1
  distinct(qName, .keep_all = T) 

# Find what chromosome has the most matches (total bp length)
match.chain.df <- 
  chain.ens.out %>% 
  group_by(qName, refName) %>% 
  summarise(MaxMatch = sum(qLength),
            refStart = min(refStart))

# Make a dataframe of maxmimum alignment summary and the best match direction for that scaffold
df.chain.raw <- 
  df.orientation %>% 
  select(-n.dir) %>% 
  right_join(match.chain.df, by = "qName")


df.chain.raw %>% filter(refName == "SCAF_11") %>% arrange(desc(MaxMatch))


# Remove scaffolds that had a cumultive matching sum < 100kb
# Format data
df.formatted <- 
  df.chain.raw %>%
  filter(MaxMatch > 50000) %>%
  dplyr::rename(SCAFF = qName, ORIENTATION = dir) 
  
old_genome_key <- read_csv(here("data_processing", "chromosemble", "outputs", "chromosome_labels.csv")) %>% rename(old_scaff = SCAFF, old_chr = CHR)
new_genome_key <- read_csv(here("data_processing", "chromosemble2", "outputs","chromosome_labels_p.csv")) %>% rename(new_scaff = scaffold, new_chr = chr)


df.formatted %>%
  group_by(refName) %>% 
  filter(MaxMatch == max(MaxMatch)) %>%
  arrange(desc(MaxMatch))

old_vs_new <-
  df.formatted %>% 
  arrange(desc(MaxMatch)) %>% 
  select(SCAFF, refName, MaxMatch) %>%
  rename(
    old_chr = SCAFF,
    new_scaff = refName
  ) %>%
  left_join(
    old_genome_key,
    by = "old_chr"
  ) %>%
  left_join(
    new_genome_key,
    by = "new_scaff"
  )  %>%
  mutate(
    old_scaff = case_when(is.na(old_scaff) ~ old_chr, TRUE ~ old_scaff)
  )%>%
  select(
    old_chr,
    new_scaff,
    new_chr,
    MaxMatch
  )

old_vs_new

chr <- c(paste0("chr", 1:11), "XY?", "Scaffold_13__1_contigs__length_49873245")
old_vs_new %>% filter(old_chr %in% chr) 
old_vs_new %>% filter(new_chr %in% chr) %>% group_by(new_chr) %>% filter(MaxMatch == max(MaxMatch)) %>% arrange(desc(MaxMatch))


old_vs_new %>% print(n = 15)
old_vs_new %>% 
  write_csv(here("data_processing", "chromosemble2", "outputs", "old_vs_new.csv"))


old_vs_new %>% filter(old_chr == "chr10")

old_vs_new %>% filter(MaxMatch > 1000000) %>% write_csv(here("data_processing", "chromosemble2", "outputs", "old_vs_new.csv"))
