# Run in analysis/anne

# Get callable sites bed
CALLABLE=../../data/ccgp_data/58-Sceloporus_callable_sites.bed
# get vcfs
VCF=../../data/ccgp_data/58-Sceloporus_complete_coords_annotated.vcf.gz
INDELS=../../data/ccgp_data/58-Sceloporus_clean_indels.vcf.gz

# Check the start of chr6 in callable sites bed file; shows that there are regions in the interval 1-24 Mbp
awk '$1 == "chr6"' $CALLABLE | head -n 1
awk '$1 == "chr6" && $2 >= 0 && $3 <= 24000000' $CALLABLE
# First pos in clean snps VCF for chr6 is 2451963
tabix $VCF chr6 | head -n 1 | awk '{print $1, $2}'
# First pos in indels VCF for chr6 is 24520191
bcftools query -f '%POS\n' $INDELS | wc -l # 5854339

# Create a new vcf with only chr6 sites from the complete coords annotated vcf
bcftools view -r chr6 $VCF -Oz -o chr6.vcf.gz
bcftools query -f '%POS\n' chr6.vcf.gz | wc -l # 4721650
# See lengths (?) of individual scaffolds
bcftools index -s $VCF
# chr6    143979556       4721650

# Count chr6 callable sites length in R
library(vcfR)
library(here)
library(tidyverse)

callbed <- read_tsv(here("data", "ccgp_data", "58-Sceloporus_callable_sites.bed"), col_names = c("locus", "startPos", "endPos"))
chr6 <- callbed %>% filter(grepl("chr6", locus)) %>% mutate(length = endPos-startPos)
chr6 %>% summarize(sum = sum(length)) # 113,360,742
