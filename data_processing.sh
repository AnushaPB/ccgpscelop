#conda env create -f ccgpscelop.yml
source activate ccgpscelop

# filter 1 (site based filter of maf and depth)
vcftools --gzvcf 58-Sceloporus_clean_snps.vcf.gz --minDP 5 --maxDP 50 --maf 0.05 --recode --recode-INFO-all --out 58-Sceloporus_maf05_minDP5_maxDP50

# visualize
vcftools --gzvcf 58-Sceloporus_maf05_minDP5_maxDP50.recode.vcf --depth --out sample_depth_info
vcftools --gzvcf 58-Sceloporus_maf05_minDP5_maxDP50.recode.vcf --missing-indv --out sample_missing_info

vcftools --gzvcf 58-Sceloporus_maf05_minDP5_maxDP50.recode.vcf --site-mean-depth --out site_depth_info
vcftools --gzvcf 58-Sceloporus_maf05_minDP5_maxDP50.recode.vcf --missing-site --out site_missing_info

# MISSING FILTER --------------------------------------------------------------------
# filter out individuals with missingness >60% (this drops one individual)
awk '$5 <= 0.60 {print $1}' sample_missing_info.imiss > keep.list
vcftools --vcf 58-Sceloporus_maf05_minDP5_maxDP50.recode.vcf --keep keep.list --recode --out 58-Sceloporus_maf05_minDP5_maxDP50_rmsamp60

# filter out sites with missingness >20%
vcftools --vcf 58-Sceloporus_maf05_minDP5_maxDP50_rmsamp60.recode.vcf --max-missing 0.8 --recode --recode-INFO-all --out 58-Sceloporus_maf05_minDP5_maxDP50_rmsamp60_mm80

# count number of SNPs 
grep -vc "^#" 58-Sceloporus_maf05_minDP5_maxDP50_rmsamp60_mm80.recode.vcf #16,287,843

# get missing data again
vcftools --vcf 58-Sceloporus_maf05_minDP5_maxDP50_rmsamp60_mm80.recode.vcf --missing-indv --out sample_depth_info_postfilter

# filter out individuals with missingness <20%
# keeps 112 out of 126
awk '$5 <= 0.20 {print $1}' sample_depth_info_postfilter.imiss > keep.list
awk '$5 > 0.20 {print $1}' sample_depth_info_postfilter.imiss > remove.list
vcftools --vcf 58-Sceloporus_maf05_minDP5_maxDP50_rmsamp60_mm80.recode.vcf --keep keep.list --recode --out 58-Sceloporus_maf05_minDP5_maxDP50_rmsamp60_mm80_rmsamp20

# rename files
mv 58-Sceloporus_maf05_minDP5_maxDP50_rmsamp60_mm80_rmsamp20.recode.vcf 58-Sceloporus_maf05_minDP5_maxDP50_rmsamp60_mm80_rmsamp20.vcf

# compress files
bgzip -c 58-Sceloporus_maf05_minDP5_maxDP50_rmsamp60_mm80_rmsamp20.vcf > 58-Sceloporus_maf05_minDP5_maxDP50_rmsamp60_mm80_rmsamp20.vcf.gz

# count number of SNPs 
zgrep -vc "^#" 58-Sceloporus_maf05_minDP5_maxDP50_rmsamp60_mm80_rmsamp20.vcf.gz #16,287,843


# LINKAGE PRUNING -------------------------------------------------------------------------
# Run linkage pruning
# window size = 50
# step size = 5
# r = 0.6

# Perform linkage pruning
# 1. add SNP IDs to VCF
# 2. perform LD pruning using plink (this gives SNPs to prune)
# 3. use prune.in file to filter the vcf
bcftools annotate --set-id '%CHROM\_%POS' 58-Sceloporus_maf05_minDP5_maxDP50_rmsamp60_mm80_rmsamp20.vcf.gz -Oz -o 58-Sceloporus_maf05_minDP5_maxDP50_rmsamp60_mm80_rmsamp20_withIDs.vcf.gz
plink2 --vcf 58-Sceloporus_maf05_minDP5_maxDP50_rmsamp60_mm80_rmsamp20_withIDs.vcf.gz --make-bed --indep-pairwise 50 5 0.6 --out 58-Sceloporus_r60 --allow-extra-chr --autosome-num 95 --const-fid --bad-freqs 
bcftools view -i 'ID=@58-Sceloporus_r60.prune.in' -O z -o 58-Sceloporus_maf05_minDP5_maxDP50_rmsamp60_mm80_rmsamp20_r60.vcf.gz 58-Sceloporus_maf05_minDP5_maxDP50_rmsamp60_mm80_rmsamp20_withIDs.vcf.gz

# check number of variants
zgrep -vc "^#" 58-Sceloporus_maf05_minDP5_maxDP50_rmsamp60_mm80_rmsamp20_r60.vcf.gz  #7,874,541/16,287,843 variants removed = 8,413,302 remaining

# CALCULATE GENETIC DISTANCE -------------------------------------------------------------------------
# using plink 1.9

# calculate distance basted on allele counts
plink --vcf 58-Sceloporus_maf05_minDP5_maxDP50_rmsamp60_mm80_rmsamp20_r60.vcf.gz --out 58-Sceloporus_plinkdist --allow-extra-chr --autosome-num 95 --distance square --const-fid

# calculate distance based on 1 - IBS (1 - proportion of shared alleles)
plink --vcf 58-Sceloporus_maf05_minDP5_maxDP50_rmsamp60_mm80_rmsamp20_r60.vcf.gz --out 58-Sceloporus_plinkdist_1ibs --allow-extra-chr --autosome-num 95 --distance square 1-ibs --const-fid


bcftools query -l CCGP/58-Sceloporus_annotated.vcf.gz > samples.txt
