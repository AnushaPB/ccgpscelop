# Satsuma alignment -------------------------------------------------------
# Modified from code provided Erik Enbody 
library(tidyverse)
library(here)

# Load chain file
chain.ens <- read_table(here("data_processing", "chromosemble", "outputs", "xcorr_first", "satsuma_summary.chained.out"), col_names = FALSE)
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

# Remove scaffolds that had a cumultive matching sum < 5kb
df.chain.filtered <- 
  arrange(df.chain.raw, refName, refStart, MaxMatch) %>% 
  mutate(refName = ifelse(MaxMatch < 5000, NA, refName))

# Select the best scaffold/chromosome matching, i.e. the one with the longest possible alignment
df.good.matches <- 
  df.chain.filtered %>% 
  group_by(qName) %>% 
  filter(MaxMatch == max(MaxMatch)) %>% 
  arrange(qName) %>% 
  filter(n() == 1)
  
# Sometimes there are two that are equal best alignments, especially for really tiny scaffolds
# Probably these can just be dropped but for now I give them chromosome NA
df.multimatches <- 
  df.chain.filtered %>% 
  group_by(qName) %>% 
  filter(MaxMatch == max(MaxMatch)) %>% arrange(qName) %>% 
  filter(n() > 1) %>% 
  mutate(refName = "NA") %>% 
  distinct(qName, .keep_all = T)

# Put the above two together
df.chain <- rbind(df.good.matches, df.multimatches)

# Format data
df.formatted <- 
  df.chain %>%
  dplyr::rename(SCAFF = qName, ORIENTATION = dir) %>% 
  # Label non-chromosomes
  mutate(refName = case_when(
    grepl("chr", refName) ~ refName,
    TRUE ~ "chrunknown"
  ))

df.formatted%>% filter(SCAFF == "Scaffold_13__1_contigs__length_49873245") %>% pull(refName)

# Get chromosomes
chromosomes <- 
  df.formatted %>% 
  filter(refName != "chrunknown") %>%
  # Transform "chromosome_X" to "chrX" and remove everything after the number
  mutate(refName = sub(".*chromosome_(\\d+).*", "chr\\1", refName)) %>%
  select(SCAFF, refName, MaxMatch)

# Get rows where chrosome is duplicated (i.e., multiple scaffolds map to the same chromosome)
# MAKE SURE TO CHECK THESE FOR MULTIPLE LARGE SCAFFOLDS ALIGNING TO THE SAME CHROSOME
dups <- chromosomes %>% ungroup() %>% count(refName) %>% filter(n > 1) %>% pull(refName)
chromosomes %>% filter(refName %in% dups) %>% arrange(refName) %>% data.frame()

# Identify the scaffold with the most matches for each chromosome
final_chromosomes <- 
  chromosomes %>% 
  group_by(refName) %>% 
  filter(MaxMatch == max(MaxMatch)) %>% 
  ungroup() 

head(final_chromosomes)

# Create output file
out <- 
  final_chromosomes %>%
  rename(CHR = refName) %>%
  mutate(CHR = factor(CHR, levels = paste0("chr", 1:nrow(final_chromosomes)))) %>%
  arrange(CHR) %>%
  select(-MaxMatch) 

# Write output
write_csv(out, here("data_processing", "chromosemble", "outputs", "chromosome_labels.csv"))
