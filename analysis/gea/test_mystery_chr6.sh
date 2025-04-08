# Get callable sites bed
CALLABLE=../../data/ccgp_data/58-Sceloporus_callable_sites.bed
# get vcf
VCF=../../data/ccgp_data/58-Sceloporus_clean_snps.vcf.gz

# Look at first positions
awk '$1 == "chr6"' $CALLABLE | head -n 1
tabix $VCF chr6 | head -n 1 | awk '{print $1, $2}'

