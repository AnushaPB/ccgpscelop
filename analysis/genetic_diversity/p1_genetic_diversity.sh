# Run in analysis/genetic_diversity
source activate ccgpscelop
BASE_PATH=../../data/ccgp_data
# !!! NOTE: running into issues with this file so currently het/window_pi was created on Rancor using Anne's 58-Sceloporus_annotated plink files
PLINK=$BASE_PATH/58-Sceloporus_complete_coords_annotated_chr
PLINK_PRUNED=$BASE_PATH/58-Sceloporus_annotated_pruned_0.6_chr

# Calculate how many SNPs in PLINK
wc -l ${PLINK}.bim #61896673

# Count how many individuals in PLINK file
wc -l ${PLINK}.fam # 257
wc -l ${PLINK_PRUNED}.fam # 257

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
# Calculate heterozygosity for only chromosomes
plink --bfile $PLINK --het --out outputs/58-Sceloporus_chr --allow-extra-chr

# Calculate heterozygosity from pruned data
plink --bfile $PLINK_PRUNED --het --out outputs/58-Sceloporus_pruned --allow-extra-chr

# Calculate heterozygosity from GEA data
plink --bfile ../gea/outputs/nonsyn --het --out outputs/nonsyn --allow-extra-chr
