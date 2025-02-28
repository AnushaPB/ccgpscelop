# Run in analysis/genetic_diversity
source activate ccgpscelop
BASE_PATH=../../data/ccgp_data
# !!! NOTE: running into issues with this file so currently het/window_pi was created on Rancor using Anne's 58-Sceloporus_annotated plink files
PLINK=$BASE_PATH/58-Sceloporus_complete_coords_annotated_chr
PLINK_PRUNED=$BASE_PATH/58-Sceloporus_complete_coords_pruned_0.6_chr

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

# Calculate heterozygosity from pruned data
plink --bfile $PLINK_PRUNED --het --out outputs/58-Sceloporus_pruned --allow-extra-chr

# PCA ---------------------------------------------------------------------------------
plink --bfile $PLINK_PRUNED --pca 3 --out outputs/58-Sceloporus_pca