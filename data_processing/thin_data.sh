conda activate ccgpscelop

PLINK=../data/ccgp_data/58-Sceloporus_complete_coords_annotated

plink --bfile $PLINK --bp-space 10000 --make-bed --allow-extra-chr --out ../data/ccgp_data/58-Sceloporus_thinned10kb

plink --bfile $PLINK --bp-space 1000 --make-bed --allow-extra-chr --out ../data/ccgp_data/58-Sceloporus_thinned1kb