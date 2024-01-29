#conda env create -f ccgpscelop.yml
source activate ccgpscelop

# filter 1
vcftools --gzvcf 58-Sceloporus_clean_snps.vcf.gz --minDP 5 --maxDP 50 --maf 0.05 --recode --recode-INFO-all --out 58-Sceloporus_maf05_minDP5_maxDP50

# visualize
vcftools --gzvcf 58-Sceloporus_maf05_minDP5_maxDP50.recode.vcf --depth --out sample_depth_info
vcftools --gzvcf 58-Sceloporus_maf05_minDP5_maxDP50.recode.vcf --missing-indv --out sample_missing_info

vcftools --gzvcf 58-Sceloporus_maf05_minDP5_maxDP50.recode.vcf --site-mean-depth --out site_depth_info
vcftools --gzvcf 58-Sceloporus_maf05_minDP5_maxDP50.recode.vcf --missing-site --out site_missing_info

# filter 2
vcftools --vcf 58-Sceloporus_maf05_minDP5_maxDP50.recode.vcf --max-missing 0.8 --recode --recode-INFO-all --out 58-Sceloporus_maf05_minDP5_maxDP50_mm80

# rename files
mv 58-Sceloporus_maf05_minDP5_maxDP50_mm80.recode.vcf 58-Sceloporus_maf05_minDP5_maxDP50_mm80.vcf

# compress files
bgzip -c 58-Sceloporus_maf05_minDP5_maxDP50_mm80.vcf > 58-Sceloporus_maf05_minDP5_maxDP50_mm80.vcf.gz
