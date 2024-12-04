source activate ccgpscelop
BASE_PATH=../../data/ccgp_data
# !!! NOTE: running into issues with this file so currently het/window_pi was created on Rancor using Anne's 58-Sceloporus_annotated plink files
PLINK=$BASE_PATH/58-Sceloporus
VCF=$BASE_PATH/58-Sceloporus_clean_snps.vcf.gz

# HETEROZYGOSITY ----------------------------------------------------------------------
# note: outputs homozygosity information
# FID: Family ID
# IID: Individual ID
# O(HOM): Observed number of homozygous genotypes
# O(HET): Observed number of heterozygous genotypes
# N(NM): Count of non-missing genotypes
# F: Inbreeding coefficient estimate
# calculate heterozygosity stats
# set const-fid to set FID (population ID) to 0; otherwise Error: Multiple instances of '_' in sample ID.
plink --bfile $PLINK --het --out outputs/58-Sceloporus --allow-extra-chr

# WINDOWED PI -------------------------------------------------------------------------
vcftools --gzvcf $VCF --window-pi 10000 --out 58-Sceloporus_10kb_windowpi
