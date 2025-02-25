# GETTING CHECKSUMS ERROR MESSAGES FOR BOTH THE FOLLOWING:
bcftools query -f '%POS\n' 58-Sceloporus_complete_coords_annotated.vcf.gz | wc -l # 7201145
bcftools query -f '%POS\n' 58-Sceloporus_complete_coords_pruned_0.6.vcf.gz | wc -l # 40705025

# Erik's call for CCGP:
# plink2 --vcf {input.vcf} --make-bed --indep-pairwise 50 5 0.6 --out {params.prefix2} --allow-extra-chr --autosome-num 95 --const-fid --bad-freqs

# Anne's call for retrieving linked SNPs:
# PLINK=../../data/ccgp_data/58-Sceloporus_complete_coords_annotated
# plink --bfile $PLINK --r2 --ld-window-kb 50 --ld-window-r2 0.6 --allow-extra-chr --autosome-num 95 --const-fid --out outputs/snp_r2

# plink2 --vcf {input.vcf} --make-bed --indep-pairwise 50 5 0.6 --out {params.prefix2} --allow-extra-chr --autosome-num 95 --const-fid --bad-freqs

mamba activate ccgpscelop
PLINK=../../data/ccgp_data/58-Sceloporus_complete_coords_annotated
# 71154824 variants and 161 samples pass filters and QC
plink --bfile $PLINK --r2 --ld-window-kb 50 --ld-window-r2 0.6 --ld-window 6 --allow-extra-chr --autosome-num 95 --const-fid --out plink_r2
plink --bfile $PLINK --r2 --ld-window-kb 50 --ld-window-r2 0.6 --ld-window 5 --allow-extra-chr --autosome-num 95 --const-fid --out plink_r2_5

PLINK=../../data/ccgp_data/58-Sceloporus_complete_coords_annotated
plink --bfile $PLINK --r2 --ld-window-kb 50 --ld-window-r2 0.6 --allow-extra-chr --autosome-num 95 --const-fid --out plink_r2

# r <- read_table(here("analysis", "anne", "plink_r2.ld"))
# nrow(r)
# [1] 9133229

r_new <- data.table::fread(here("analysis", "anne", "plink_r2_5.ld"))
nrow(r_new)
# [1] 7748715

r_new %>% filter((SNP_A == "chr1_3757363_G_A" & SNP_B == "chr1_3757379_T_C") | (SNP_A == "chr1_3757379_T_C" & SNP_B == "chr1_3757363_G_A"))

### =========== TRY PLINK2 AGAIN

# New LD pruning line run by Erik:
plink2 --vcf {input.vcf} --make-bed --indep-pairwise 50kb 1 0.6 --out {params.prefix2} --allow-extra-chr --autosome-num 95 --const-fid --bad-freqs

# Make plink2's bfiles
plink2 --vcf 58-Sceloporus_complete_coords_annotated.vcf.gz --make-bed --allow-extra-chr --autosome-num 95 --const-fid --bad-freqs --out 58-Sceloporus_complete_coords_annotated

PLINK=../../data/ccgp_data/58-Sceloporus_complete_coords_annotated
plink2 --bfile $PLINK --r2-unphased inter-chr --ld-window-kb 50 --ld-window-r2 0.6 --allow-extra-chr --autosome-num 95 --const-fid --bad-freqs --out plink2_r2

PLINK=../../data/ccgp_data/58-Sceloporus_complete_coords_annotated
# When a limited window report is requested, every pair of variants with at least (10-1) variants between them, or more than 1000 kilobases apart, is ignored. You can change the first threshold with --ld-window, and the second threshold with --ld-window-kb.
plink --bfile $PLINK --r2 --ld-window 2 --ld-window-kb 50 --ld-window-r2 0.6 --allow-extra-chr --autosome-num 95 --const-fid --out plink_r2

### Because we can't calculate backwards correlations with plink2, let's generate an LD-pruning report using plink1.9 and see how that compares
