# Run in analysis/admixture
source activate ccgpscelop
BASE_PATH=../../data/ccgp_data
PLINK_PRUNED=$BASE_PATH/58-Sceloporus_annotated_pruned_0.6_chr

plink --bfile $PLINK_PRUNED --pca 3 --out outputs/58-Sceloporus_pca