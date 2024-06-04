# run in data/processed_data
#conda env create -f ccgpscelop.yml
source activate ccgpscelop
RAW_DATA=../raw_data/QC/58-Sceloporus_filteredQC.vcf.gz

# REMOVE: Scelocci_CHI1382_DAW5-46-21 from file (uncertain  provenance - QC plots suggest location is incorrect)
vcftools --remove-indv Scelocci_CHI1382_DAW5-46-21 --gzvcf $RAW_DATA --recode --out 58-Sceloporus_filteredQC_samplerm

# Reassign raw data object
RAW_DATA=58-Sceloporus_filteredQC_samplerm.recode.vcf

# calculate missingness and depth data from raw data
vcftools --gzvcf $RAW_DATA --depth --out sample_depth_info_raw
vcftools --gzvcf $RAW_DATA --missing-indv --out sample_missing_info_raw

vcftools --gzvcf $RAW_DATA --site-mean-depth --out site_depth_info_raw
vcftools --gzvcf $RAW_DATA --missing-site --out site_missing_info_raw

# filter 1 (site based filter of maf and depth)
vcftools --gzvcf $RAW_DATA --minDP 5 --maxDP 50 --maf 0.05 --recode --recode-INFO-all --out 58-Sceloporus_maf05_minDP5_maxDP50

# calculate missingness and depth data
vcftools --gzvcf 58-Sceloporus_maf05_minDP5_maxDP50.recode.vcf --depth --out sample_depth_info
vcftools --gzvcf 58-Sceloporus_maf05_minDP5_maxDP50.recode.vcf --missing-indv --out sample_missing_info

vcftools --gzvcf 58-Sceloporus_maf05_minDP5_maxDP50.recode.vcf --site-mean-depth --out site_depth_info
vcftools --gzvcf 58-Sceloporus_maf05_minDP5_maxDP50.recode.vcf --missing-site --out site_missing_info

# MISSING FILTER --------------------------------------------------------------------
# filter out individuals with missingness >60% (this drops 2 individuals)
awk '$5 <= 0.60 {print $1}' sample_missing_info.imiss > keep.list
awk '$5 > 0.60 {print $1}' sample_missing_info.imiss > remove.list
vcftools --vcf 58-Sceloporus_maf05_minDP5_maxDP50.recode.vcf --keep keep.list --recode --out 58-Sceloporus_maf05_minDP5_maxDP50_rmsamp60

# filter out sites with missingness >20%
vcftools --vcf 58-Sceloporus_maf05_minDP5_maxDP50_rmsamp60.recode.vcf --max-missing 0.8 --recode --recode-INFO-all --out 58-Sceloporus_maf05_minDP5_maxDP50_rmsamp60_mm80

# count number of SNPs 
grep -vc "^#" 58-Sceloporus_maf05_minDP5_maxDP50_rmsamp60_mm80.recode.vcf #14877047

# get missing data again
vcftools --vcf 58-Sceloporus_maf05_minDP5_maxDP50_rmsamp60_mm80.recode.vcf --missing-indv --out sample_depth_info_postfilter

# filter out individuals with missingness >40% (this drops 2 individuals)
# Note: the max missingness is actually 0.32 with this cut-off
awk '$5 <= 0.40 {print $1}' sample_depth_info_postfilter.imiss > keep_postfilter.list
awk '$5 > 0.40 {print $1}' sample_depth_info_postfilter.imiss > remove_postfilter.list
vcftools --vcf 58-Sceloporus_maf05_minDP5_maxDP50_rmsamp60_mm80.recode.vcf --keep keep_postfilter.list --recode --out 58-Sceloporus_maf05_minDP5_maxDP50_rmsamp60_mm80_rmsamp40

# Set prefix
PREFIX=58-Sceloporus_maf05_minDP5_maxDP50_rmsamp60_mm80_rmsamp40

# rename files
cp $PREFIX.recode.vcf $PREFIX.vcf

# compress files
bgzip -c $PREFIX.vcf > $PREFIX.vcf.gz

# count number of SNPs 
zgrep -vc "^#" $PREFIX.vcf.gz 

# LINKAGE PRUNING -------------------------------------------------------------------------
# Run linkage pruning
# window size = 50
# step size = 5
# r = 0.6

# Perform linkage pruning
# 1. add SNP IDs to VCF
# 2. perform LD pruning using plink (this gives SNPs to prune)
# 3. use prune.in file to filter the vcf
bcftools annotate --set-id '%CHROM\_%POS' $PREFIX.vcf.gz -Oz -o ${PREFIX}_withIDs.vcf.gz
plink2 --vcf ${PREFIX}_withIDs.vcf.gz --make-bed --indep-pairwise 50 5 0.6 --out 58-Sceloporus_r60 --allow-extra-chr --autosome-num 95 --const-fid --bad-freqs 
bcftools view -i 'ID=@58-Sceloporus_r60.prune.in' -O z -o ${PREFIX}_r60.vcf.gz ${PREFIX}_withIDs.vcf.gz

# check number of variants
zgrep -vc "^#" ${PREFIX}_r60.vcf.gz #7,643,887

# create plink files (used for FEEMS and ADMIXTURE)
plink2 --vcf ${PREFIX}_r60.vcf.gz --make-bed --out ${PREFIX}_r60 --allow-extra-chr --autosome-num 95 --const-fid

# CALCULATE GENETIC DISTANCE -------------------------------------------------------------------------
# using plink 1.9

# calculate distance basted on allele counts
plink --vcf ${PREFIX}_r60.vcf.gz --out 58-Sceloporus_plinkdist --allow-extra-chr --autosome-num 95 --distance square --const-fid

# calculate distance based on 1 - IBS (1 - proportion of shared alleles)
plink --vcf ${PREFIX}_r60.vcf.gz --out 58-Sceloporus_plinkdist_1ibs --allow-extra-chr --autosome-num 95 --distance square 1-ibs --const-fid

# get samples
bcftools query -l CCGP/58-Sceloporus_annotated.vcf.gz > samples.txt

# FILTER ONE CONTIG -----------------------------------------------------------------------------------
# get first contig name: JALMGF010000001.1
# zgrep "^##contig=" ${PREFIX}_r60.vcf.gz | head -n 1 | awk -F'[=,]' '{print $3}'

# index file
tabix -p vcf ${PREFIX}_r60.vcf.gz

# filter out that contig
bcftools view -r JALMGF010000001.1 ${PREFIX}_r60.vcf.gz > 58-Sceloporus_JALMGF010000001.1.vcf.gz

# COUNT NUMBER OF SITES ------------------------------------
# count snps
zgrep -vc "^#" $RAW_DATA > 58-Sceloporus_allsnps_nsites.txt

# count snps postfilters
zgrep -vc "^#" $PREFIX.vcf.gz > 58-Sceloporus_postfiltersnps_nsites.txt

# count number of callable sites 
python ../../data_processing/count_callable.py

# COPY TO RANCOR --------------------------------------------
#rclone copy 58-Sceloporus_maf05_minDP5_maxDP50_rmsamp60_mm80_rmsamp40.vcf.gz rancor:
#rclone copy 58-Sceloporus_maf05_minDP5_maxDP50_rmsamp60_mm80_rmsamp40_r60.vcf.gz rancor:


# ALLELE COUNTS ---------------------------------------------
# Create maf < 0.01 data with just the samples that passed the filters
RAW_DATA=../raw_data/QC/58-Sceloporus_filteredQC.vcf.gz
PREFIX2=58-Sceloporus_maf01_minDP10_maxDP50
vcftools --gzvcf $RAW_DATA --minDP 10 --maxDP 50 --keep keep_postfilter.list --recode --recode-INFO-all --out $PREFIX2

# Get allele counts
bcftools +fill-tags $PREFIX2.recode.vcf -Ov -o filled.vcf -- -t AC,AN
bcftools query -f'[%CHROM\t%POS\t%REF\t%ALT\t%AC\t%AN\n]' filled.vcf | uniq > allele_counts.txt

# Filter allele counts
awk '$5 < 4' allele_counts.txt > filtered_allele_counts.txt

# Get length
wc -l filtered_allele_counts.txt 

# Get positions
awk '{print $1"\t"$2"\t"$3"\t"$4}' filtered_allele_counts.txt > positions.list

# Use bcftools to view the variants of interest
# from: https://www.biostars.org/p/170965/
awk 'BEGIN{print "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"} {print $1"\t"$2"\t.\t"$3"\t"$4"\t.\t.\t."}' positions.list > positions.vcf
bgzip -c positions.vcf > positions.vcf.gz
tabix -p vcf positions.vcf.gz
bgzip -c $PREFIX2.recode.vcf > $PREFIX2.vcf.gz
tabix -p vcf $PREFIX2.vcf.gz
bcftools isec -n=2 -w1 -O v -o rare_alleles.vcf $PREFIX2.vcf.gz positions.vcf.gz

# Check number of sites
grep -vc "^#" rare_alleles.vcf 