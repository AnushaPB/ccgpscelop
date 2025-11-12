# Run in analysis/ibdibe
source activate ccgpscelop
BASE_PATH=../../data/ccgp_data
# !!! NOTE: running into issues with this file so currently het/window_pi was created on Rancor using Anne's 58-Sceloporus_annotated plink files
PLINK=$BASE_PATH/58-Sceloporus_annotated_pruned_0.6_chr

plink --bfile $PLINK --allow-extra-chr --distance square 1-ibs --const-fid --out outputs/58-Sceloporus_annotated_pruned_0.6_chr

plink --bfile ../gea/outputs/nonsyn --allow-extra-chr --distance square 1-ibs --const-fid --out outputs/nonsyn

plink --bfile ../gea/outputs/genes --allow-extra-chr --distance square 1-ibs --const-fid --out outputs/genes

plink --bfile ../gea/outputs/gea --allow-extra-chr --distance square 1-ibs --const-fid --out outputs/gea