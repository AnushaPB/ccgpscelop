PLINK=../data/ccgp_data/58-Sceloporus_annotated_complete_coords

plink --bfile $PLINK --bp-space 10000 --make-bed --allow-extra-chr --out ../data/ccgp_data/58-Sceloporus_thinned