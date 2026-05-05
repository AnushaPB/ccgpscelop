library(pafr)
# library(tidyverse)
# # https://cran.r-project.org/web/packages/pafr/vignettes/Introduction_to_pafr.html
# Column	Name	Data Type	Description
# 1	qname	string	Query sequence name
# 2	qlen	int	Query sequence length
# 3	qstart	int	Query start coordinate (0-based)
# 4	qend	int	Query end coordinate (0-based)
# 5	strand	char	‘+’ if query/target on the same strand; ‘-’ if opposite
# 6	tname	string	Target sequence name
# 7	tlen	int	Target sequence length
# 8	tstart	int	Target start coordinate on the original strand
# 9	tend	int	Target end coordinate on the original strand
# 10	nmatch	int	Number of matching bases in the mapping
# 11	alen	int	Number of bases, including gaps, in the mapping
# 12	mapq	int	Mapping quality (0-255, with 255 if missing)

paf <- read_table("map_GCA_029215755.1_aSpeHam1.0.hap2_genomic_to_GCA_023333645.1_rSceOcc1.0.p_genomic_backup.paf", col_names = FALSE)
colnames <- c("qname", "qlen", "qstart", "qend", "strand", "tname", "tlen", "tstart", "tend", "nmatch", "alen", "mapq")
names(paf) <- colnames

# Drop no name columns
paf <- paf %>% select(all_of(colnames))

# Check total sequence length
paf %>% distinct(tname, tlen) %>% pull(tlen) %>% sum() / 1e9 # 2.8 Gb

byscaff <-
  paf %>% 
  filter(mapq > 20) %>%
  # You can check results by adding the following filter and confirming that JARDYJ010000005.1 is what is returned
  #filter(tname == "JALMGF010000011.1") %>% 
  group_by(qname, tname) %>%
  summarize(alen = sum(alen), nmatch = sum(nmatch), tlen = unique(tlen)) %>%
  mutate(perc = alen/tlen) %>%
  arrange(desc(perc)) %>%
  group_by(tname) %>%
  filter(perc == max(perc)) 

# PROBLEM: percent > 1 because of multiple alignments to same region
# Example, the tlen is < than the summed alen (note that tlen is short in general):
paf %>% 
  filter(mapq > 20) %>%
  filter(tname == "JALMGF010000142.1") %>%
  select(qname, qstart, qend, tname, tlen, alen, nmatch) 

paf %>% 
  filter(tname == "JALMGF010000142.1") %>%
  select(mapq, qname, qstart, qend, tname, tlen, alen, nmatch) %>%
  mutate(alen/tlen)

byscaff %>% filter(tname == "JALMGF010000142.1")

# Filtering out alignments that are 90% of the length of the reference
byscaff90 <- byscaff %>% filter(perc > 0.10)

# I trust the longer tlens
byscaff90 %>% arrange(desc(tlen))

# Calc total number of scaffolds affected
byscaff90 %>% distinct(tname) %>% nrow()

# Sum total length of affected scaffolds
byscaff90 %>% distinct(tname, tlen) %>% pull(tlen) %>% sum() / 1e9 # 0.97 Gb

# Write scaffolds to file
byscaff90 %>% 
  arrange(desc(tlen)) %>%
  distinct(tname, tlen) %>% 
  rename(scaffold = tname, length = tlen) %>%
  write_tsv("occidentalis_hammondii_overlapping_scaffolds.txt")

library(tidyverse)
overlap <- read_table("occidentalis_hammondii_overlapping_scaffolds.txt")
zgs <- read_csv("zero_genotype_scaffolds.csv")
lgs <- read_csv("low_genotype_scaffolds.csv")
cgs <- bind_rows(zgs, lgs)
write_tsv(cgs, "low_or_no_genotype_scaffolds.txt")

# Get all scaffolds with low/no genotypes that are not in the contaminated set
cgs %>% filter(!(scaffold %in% overlap$scaffold)) # 81 
# Get all contaminated scaffolds that are not in the low/no genotypes
overlap %>% filter(!(scaffold %in% cgs$scaffold)) # None!

overlap %>% left_join(cgs)
mean(overlap$scaffold %in% cgs$scaffold) 
mean(cgs$scaffold %in% overlap$scaffold)

# NOT SURE ABOUT STUFF BELOW:

# CHECKING:
paf %>% 
  filter(mapq > 30) %>%
  # You can check results by adding the following filter and confirming that JARDYJ010000005.1 is what is returned
  filter(tname == "JALMGF010000011.1") %>% 
  arrange(desc(alen))%>%
   mutate(perc = alen/tlen) %>%
  select(qname, tname, tlen, alen, nmatch, perc, qstart, qend) 

# Und vs Occ
paf <- read_table("map_GCF_019175285.1_SceUnd_v1.1_genomic_to_jordan-uni4378-mb-hirise-65qvq__12-23-2023__final_assembly_relabelled.paf", col_names = FALSE)
colnames <- c("qname", "qlen", "qstart", "qend", "strand", "tname", "tlen", "tstart", "tend", "nmatch", "alen", "mapq")
names(paf) <- colnames

# Drop no name columns
paf <- paf %>% select(all_of(colnames))

# Check total sequence length
paf %>% distinct(tname, tlen) %>% pull(tlen) %>% sum() / 1e9 # 1.8 Gb
paf %>% distinct(qname, qlen) %>% pull(qlen) %>% sum() / 1e9 # 1.9 Gb

paf %>% filter(mapq > 30, nmatch > 100000, tname == "chr6") %>% distinct(qname)
paf %>% filter(mapq > 30, nmatch > 100000, qname == "NC_056527.1") %>% distinct(tname)

chr6_scaff13_occ_len <- paf %>% filter(mapq > 30, nmatch > 100000, qname == "NC_056527.1") %>% distinct(tname, tlen) %>% pull(tlen) %>% sum() # 193852802
chr6_und_len <- paf %>% filter(qname == "NC_056527.1") %>% distinct(qname, qlen) %>% pull(qlen) #170077623
chr6_scaff13_occ_len - chr6_und_len
