### Linking GEA outliers to pruned sites
# Run within analysis/anne directory
mamba activate ccgpscelop

# Get correlation values
# Specify window size and r2 value the same as our initial LD-pruning prior to running GEA
plink --vcf ../../data/ccgp_data/58-Sceloporus_annotated.vcf.gz --r2 --ld-window-kb 50 --ld-window-r2 0.6 --allow-extra-chr --autosome-num 95 --const-fid --out outputs/58-Sceloporus_snp_correlations

# Concatenate RDA results into a single file
awk 'FNR>1' ../../outputs/RDA/*/*Zscores.csv > outputs/58-Sceloporus_Zscores.csv
awk 'FNR>1' ../../outputs/RDA/*/*rdadapt.csv > outputs/58-Sceloporus_rdadapt.csv

# Retrieve SNPs correlated with RDA outliers
# rsync -avP hgdownload.soe.ucsc.edu::ccgp/CCGP-module/58-Sceloporus_chr/CCGP/58-Sceloporus_0.6.prune.out ./outputs/

# args are: prunedout_file / corr_file / gea_path / sig / output_format / output_path
Rscript GEA_outlier_positions.R "./outputs/58-Sceloporus_0.6.prune.out" "./outputs/58-Sceloporus_snp_correlations.ld" "../../outputs/RDA/*/" 0.01 "together" "./outputs/"