VCF=../data/ccgp_data/58-Sceloporus_clean_snps.vcf.gz

# Convert VCF to PLINK
plink --vcf $VCF --make-bed --out ../data/ccgp_data/58-Sceloporus --allow-extra-chr --const-fid