conda activate ccgpscelop

VCF=../data/ccgp_data/58-Sceloporus_complete_coords_annotated.vcf.gz
PRUNED=../data/ccgp_data/58-Sceloporus_complete_coords_pruned_0.6.vcf.gz

# Convert VCF to PLINK
plink --vcf $VCF --make-bed --out ../data/ccgp_data/58-Sceloporus_complete_coords_annotated --allow-extra-chr --const-fid > vcf_to_plink.out 2> vcf_to_plink.err

plink --vcf $PRUNED --make-bed --out ../data/ccgp_data/58-Sceloporus_complete_coords_pruned_0.6 --allow-extra-chr --const-fid > vcf_to_plink.out 2> vcf_to_plink.err

#plink --vcf 58-Sceloporus_complete_coords_annotated.vcf.gz --make-bed --out 58-Sceloporus_complete_coords_annotated --allow-extra-chr --const-fid

#plink --vcf 58-Sceloporus_complete_coords_pruned_0.6.vcf.gz --make-bed --out 58-Sceloporus_complete_coords_pruned_0.6 --allow-extra-chr --const-fid