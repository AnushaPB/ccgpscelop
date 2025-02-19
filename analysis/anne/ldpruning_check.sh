### Checking sig SNPs

outpath <- here("analysis", "gea", "outputs")rdasig <- read_csv(here(outpath, "RDA_bio1_ndvi", "58-Sceloporus_RDA_outliers_full_rdadapt.csv")) cortest <- read_csv(here(outpath, "RDA_bio1_ndvi", "58-Sceloporus_RDA_cortest_full.csv"))
rdasig <- read_csv(here(outpath, "RDA_bio1_ndvi", "58-Sceloporus_RDA_outliers_full_rdadapt.csv")) 
rdasig01 <- rdasig %>% filter(q.values < 0.01)
cortest_p <- cortest %>% filter(outlier_method == "p")
all(rdasig01$locus %in% cortest_p$snp) # TRUE
nrow(rdasig01)# [1] 1328466
nrow(cortest_p)# [1] 2656932

### Checking LD pruning

r <- read_table(here("analysis", "gea", "outputs", "snp_r2.ld"))
r %>% filter((SNP_A == "chr1_3757363_G_A" & SNP_B == "chr1_3757379_T_C") | (SNP_A == "chr1_3757379_T_C" & SNP_B == "chr1_3757363_G_A"))
rdasig %>% filter(locus == "chr1_3757363_G_A" | locus == "chr1_3757379_T_C")

prunein <- read_table(here("data", "ccgp_data", "58-Sceloporus_0.6.prune.in"), col_names = "snp")
nrow(prunein)
# [1] 49306298
pruneout <- read_table(here("data", "ccgp_data", "58-Sceloporus_0.6.prune.out"), col_names = "snp")
nrow(pruneout)
# [1] 19670579
all(prunein$snp %in% pruneout$snp)
# [1] FALSE
prunein %>% filter(snp == "chr1_3757363_G_A" | snp == "chr1_3757379_T_C") # both are there
pruneout %>% filter(snp == "chr1_3757363_G_A" | snp == "chr1_3757379_T_C") # neither is there

### ===========

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