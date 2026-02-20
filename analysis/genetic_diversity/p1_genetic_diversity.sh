# Run in analysis/genetic_diversity
source activate ccgpscelop
BASE_PATH=../../data/ccgp_data
# !!! NOTE: running into issues with this file so currently het/window_pi was created on Rancor using Anne's 58-Sceloporus_annotated plink files
PLINK=$BASE_PATH/58-Sceloporus_complete_coords_annotated_chr
PLINK_PRUNED=$BASE_PATH/58-Sceloporus_annotated_pruned_0.6_chr

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


# DOSAGE FILES ------------------------------------------------------------------------
plink --bfile ../../data/ccgp_data/58-Sceloporus_pruned_0.6_thinned_10kb_chr --recode A --allow-extra-chr --out outputs/thinned # creates thinned.raw

# PCA ---------------------------------------------------------------------------------
plink --bfile $PLINK_PRUNED --pca 3 --out outputs/58-Sceloporus_pca

# Run PCA for each chromosome
for i in {1..11}
do
    plink --bfile $PLINK_PRUNED --pca 3 --out outputs/58-Sceloporus_pca_chr$i --allow-extra-chr --chr $i
done

