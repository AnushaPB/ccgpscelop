# Run in analysis/gea
source activate ccgpscelop

# Get plink 
PLINK=../../data/ccgp_data/58-Sceloporus_complete_coords_annotated

# Calculate correlations
# Erik;s code: plink2 --vcf {input.vcf} --make-bed --indep-pairwise 50kb 0.6 --out {params.prefix2} --allow-extra-chr --autosome-num 95 --const-fid --bad-freqs
plink --bfile $PLINK --r2 --ld-window-kb 50 --ld-window-r2 0.6 --allow-extra-chr --autosome-num 95 --const-fid --out outputs/snp_r2
