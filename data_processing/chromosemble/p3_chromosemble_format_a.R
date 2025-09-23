# Satsuma alignment -------------------------------------------------------
# Modified from code provided Erik Enbody 
library(tidyverse)
library(here)

# Load chain file
chain.ens <- read_table(here("data_processing", "chromosemble", "outputs", "outputs_a", "xcorr_first", "satsuma_summary.chained.out"), col_names = FALSE)
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

# Remove scaffolds that had a cumultive matching sum < 100kb
df.formatted <- 
  df.chain.raw %>%
  filter(MaxMatch > 100000) %>%
  dplyr::rename(SCAFF = qName, ORIENTATION = dir) %>% 
  # Label non-chromosomes
  mutate(refName = case_when(
    grepl("chr", refName) ~ refName,
    TRUE ~ "chrunknown"
  ))

# Get chromosomes
chromosomes <- 
  df.formatted %>% 
  filter(refName != "chrunknown") %>%
  # Transform "chromosome_X" to "chrX" and remove everything after the number
  mutate(refName = sub(".*chromosome_(\\d+).*", "chr\\1", refName)) %>%
  select(SCAFF, refName, MaxMatch)

# Make key for chromosomes (for reference to S. undulatus)
chromosomes_key <- 
  df.formatted %>% 
  filter(refName != "chrunknown") %>%
  # Transform "chromosome_X" to "chrX" and remove everything after the number
  mutate(chr = sub(".*chromosome_(\\d+).*", "chr\\1", refName)) %>%
  mutate(nc = sub(".*(NC_\\d+\\.\\d+).*", "\\1", refName)) %>%
  ungroup() %>%
  select(nc, chr, refName) %>%
  distinct()
filter(chromosomes_key, chr == "chr6")

# Get rows where chrosome is duplicated (i.e., multiple scaffolds map to the same chromosome)
# MAKE SURE TO CHECK THESE FOR MULTIPLE LARGE SCAFFOLDS ALIGNING TO THE SAME CHROSOME
dups <- chromosomes %>% ungroup() %>% count(refName) %>% filter(n > 1) %>% pull(refName)
chromosomes %>% filter(refName %in% dups) %>% arrange(refName) %>% data.frame()

# Identify the scaffold with > 1Mb matches for each chromosome
final_chromosomes <- 
  chromosomes %>% 
  filter(MaxMatch > 1000000) %>%
  arrange(desc(MaxMatch)) 

final_chromosomes

# Create output file
out <- 
  final_chromosomes %>%
  rename(CHR = refName) %>%
  select(-MaxMatch) %>%
  rename(
    scaffold = SCAFF,
    und_chr = CHR
  ) %>%
  mutate(chr =
    case_when(
      und_chr == "chr6" & scaffold == "SCAF_7" ~ "chr7",
      scaffold == "SCAF_10" ~ "chr10",
      TRUE ~ und_chr
    )
  ) %>%
  select(-und_chr) 

out

# Write output
write_csv(out, here("data_processing", "chromosemble2", "outputs", "chromosome_labels_a.csv"))
