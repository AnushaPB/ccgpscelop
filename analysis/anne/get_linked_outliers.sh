### Linking GEA outliers to pruned sites
# Run within analysis/anne directory
mamba activate ccgpscelop

# Retrieve scaffolds in which RDA was run, ignoring header
awk -F "\"*,\"*" '{print $5}' ../../outputs/RDA/58-Sceloporus_RDA_outliers_full_rdadapt.csv | uniq | tail -n +2 > RDA_scaffolds

# Split vcf into two subsetted vcfs
vcftools --gzvcf ../../data/ccgp_data/58-Sceloporus_complete_coords_annotated.vcf.gz --chr Scaffold_1__1_contigs__length_4124712 --chr chr2 --chr chr3 --chr chr4 --chr Scaffold_7__1_contigs__length_3629712 --chr chr1 --chr chr5 --chr chr6 --chr chr7 --chr Scaffold_13__1_contigs__length_49873245 --recode --out 58-Sceloporus_complete_coords_annotated_subset1
vcftools --gzvcf ../../data/ccgp_data/58-Sceloporus_complete_coords_annotated.vcf.gz --chr chr8 --chr chr9 --chr chr11 --chr Scaffold_39__1_contigs__length_9170529 --chr Scaffold_56__1_contigs__length_6187996 --chr Scaffold_67__1_contigs__length_5618279 --chr Scaffold_83__1_contigs__length_4406559 --chr chr10 --chr Scaffold_105__1_contigs__length_3306775 --chr Scaffold_138__1_contigs__length_2391489  --recode --out 58-Sceloporus_complete_coords_annotated_subset2

# Verify that all SNPs are contained in both vcfs
# bcftools query -f 

# Get correlation values for all variants
# Specify window size and r2 value the same as our initial LD-pruning prior to running GEA
# Split vcf into two for plink calculations because keeping as one will run into memory issues
plink --vcf 58-Sceloporus_complete_coords_annotated_subset1.recode.vcf --r2 --ld-window-kb 50 --ld-window-r2 0.6 --allow-extra-chr --autosome-num 95 --const-fid --out 58-Sceloporus_snp_correlations_subset1
plink --vcf 58-Sceloporus_complete_coords_annotated_subset2.recode.vcf --r2 --ld-window-kb 50 --ld-window-r2 0.6 --allow-extra-chr --autosome-num 95 --const-fid --out 58-Sceloporus_snp_correlations_subset2

# MEMORY ISSUES WITH PLINK BELOW
# plink --vcf ../../data/ccgp_data/58-Sceloporus_complete_coords_annotated.vcf.gz --r2 --ld-window-kb 50 --ld-window-r2 0.6 --allow-extra-chr --autosome-num 95 --const-fid --out outputs/58-Sceloporus_snp_correlations

# Concatenate RDA results into a single file; this has already been done within ccgp/algatr/58-Sceloporus/RDA
# awk 'FNR>1' ../../outputs/RDA/*/*Zscores.csv > outputs/58-Sceloporus_Zscores.csv
# awk 'FNR>1' ../../outputs/RDA/*/*rdadapt.csv > outputs/58-Sceloporus_rdadapt.csv

# Retrieve SNPs correlated with RDA outliers
# rsync -avP hgdownload.soe.ucsc.edu::ccgp/CCGP-module/58-Sceloporus_chr/CCGP/58-Sceloporus_0.6.prune.out ./outputs/

# args are: prunedout_file / corr_file / gea_path / sig / output_format / output_path
Rscript GEA_outlier_positions.R "./outputs/58-Sceloporus_0.6.prune.out" "./outputs/58-Sceloporus_snp_correlations.ld" "../../outputs/RDA/*/" 0.01 "together" "./outputs/"