#conda env create -f ccgpscelop.yml
source activate ccgpscelop

vcftools --vcf 58-Sceloporus_snpsonly_rmsamp.vcf --minDP 5 --maxDP 50 --max-missing 0.8 --maf 0.05 --recode --recode-INFO-all --out 58-Sceloporus_strictfilter

# rename files
mv 58-Sceloporus_strictfilter.recode.vcf 58-Sceloporus_strictfilter.vcf

# compress files
bgzip -c 58-Sceloporus_strictfilter.vcf > 58-Sceloporus_strictfilter.vcf.gz
