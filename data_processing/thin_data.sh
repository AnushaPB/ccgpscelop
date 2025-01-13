conda activate ccgpscelop

PLINK=../data/ccgp_data/58-Sceloporus_complete_coords_pruned_0.6.vcf.gz

plink --bfile $PLINK --bp-space 10000 --make-bed --allow-extra-chr --out ../data/ccgp_data/58-Sceloporus_pruned_0.6_thinned_10kb