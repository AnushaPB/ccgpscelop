conda activate ccgpscelop

PLINK=../data/ccgp_data/58-Sceloporus_annotated_pruned_0.6_chr

# Make plink fil
plink --bfile $PLINK --bp-space 10000 --make-bed --allow-extra-chr --out ../data/ccgp_data/58-Sceloporus_pruned_0.6_thinned_10kb_chr