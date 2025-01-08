VCF=../data/ccgp_data/58-Sceloporus_complete_coords_annotated.vcf.gz

# Convert VCF to PLINK
plink --vcf $VCF --make-bed --out ../data/ccgp_data/58-Sceloporus_complete_coords_annotated --allow-extra-chr --const-fid