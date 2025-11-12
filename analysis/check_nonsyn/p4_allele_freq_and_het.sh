source activate ccgpscelop
BASE_PATH=../../data/ccgp_data
PLINK=$BASE_PATH/58-Sceloporus_complete_coords_annotated_chr

# Calculate frequencies of each allele
plink --bfile $PLINK --freq --out outputs/58-Sceloporus --allow-extra-chr
plink --bfile $PLINK --extract range outputs/all_nonsynonymous.bed --freq --out outputs/all_nonsynonymous --allow-extra-chr
plink --bfile $PLINK --extract range outputs/all_synonymous.bed --freq --out outputs/all_synonymous --allow-extra-chr

# Calculate observed heterozygosity
plink --bfile $PLINK --het --out outputs/58-Sceloporus --allow-extra-chr
plink --bfile $PLINK --extract range outputs/all_nonsynonymous.bed --het --out outputs/all_nonsynonymous --allow-extra-chr
plink --bfile $PLINK --extract range outputs/all_synonymous.bed --het --out outputs/all_synonymous --allow-extra-chr