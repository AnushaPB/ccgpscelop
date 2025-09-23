# Run in analysis/ibdibe
source activate ccgpscelop
BASE_PATH=../../data/ccgp_data
# !!! NOTE: running into issues with this file so currently het/window_pi was created on Rancor using Anne's 58-Sceloporus_annotated plink files
PLINK=$BASE_PATH/58-Sceloporus_annotated_pruned_0.6_chr
GEA_PLINK=../gea/outputs/genes

plink --bfile $PLINK --allow-extra-chr --autosome-num 95 --distance square 1-ibs --const-fid --out outputs/58-Sceloporus_annotated_pruned_0.6_chr

plink --bfile $GEA_PLINK --allow-extra-chr --autosome-num 95 --distance square 1-ibs --const-fid --out outputs/genes
