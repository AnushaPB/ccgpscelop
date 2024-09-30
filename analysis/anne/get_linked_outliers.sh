### Linking GEA outliers to pruned sites
# Run within analysis/anne directory

# Get correlation values
# Specify window size and r2 value the same as our initial LD-pruning prior to running GEA
plink --vcf ../../data/58-Sceloporus_maf05_minDP5_maxDP50_rmsamp60_mm80_rmsamp40.vcf.gz --r2 --ld-window-kb 50 --ld-window-r2 0.6 --allow-extra-chr --autosome-num 95 --const-fid --out outputs/58-Sceloporus_snp_correlations

# Concatenate RDA results into a single file
awk 'FNR>1' ../../gea/results/*/algatr/subsets/RDA/*/*Zscores.csv > outputs/58-Sceloporus_Zscores.csv
awk 'FNR>1' ../../gea/results/*/algatr/subsets/RDA/*/*rdadapt.csv > outputs/58-Sceloporus_rdadapt.csv

# Retrieve SNPs correlated with RDA outliers
# args are: prunedout_file / corr_file / gea_path / sig / output_format / output_path
Rscript GEA_outlier_positions.R "../../58-Sceloporus.prune.out" "outputs/58-Sceloporus_snp_correlations.ld" "../gea/outputs" 0.01 "together" "outputs/"