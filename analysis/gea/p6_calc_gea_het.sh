# Run in analysis/gea
source activate ccgpscelop

# Get plink 
PLINK=../../data/ccgp_data/58-Sceloporus_complete_coords_annotated

# Subset plink file with RDA snps (created using p1_process_rda.R)
plink --bfile $PLINK \
      --extract outputs/rda_ids.txt \
      --make-bed \
      --out outputs/gea --allow-extra-chr

# Subset plink file with RDA genes (created using p3_intersect_genes.R)
plink --bfile outputs/gea \
      --extract outputs/gea_gene_ids.txt \
      --make-bed \
      --out outputs/genes --allow-extra-chr

# Subset plink file without RDA snps
plink --bfile $PLINK \
      --exclude outputs/rda_ids.txt \
      --make-bed \
      --out outputs/nogea --allow-extra-chr

# Subset nogea plink file to the same number of snps as genes plink file
plink --bfile outputs/genes --write-snplist --out outputs/genes_snplist
N=$(wc -l < outputs/genes_snplist.snplist)
plink --bfile outputs/nogea --thin-count $N --out outputs/nogea_thinned --make-bed --allow-extra-chr

# Calculate heterozygosity stats
# note: outputs homozygosity information
# FID: Family ID
# IID: Individual ID
# O(HOM): Observed number of homozygous genotypes
# O(HET): Observed number of heterozygous genotypes
# N(NM): Count of non-missing genotypes
# F: Inbreeding coefficient estimate
# calculate heterozygosity stats
# set const-fid to set FID (population ID) to 0; otherwise Error: Multiple instances of '_' in sample ID.
plink --bfile outputs/gea --het --out outputs/gea --allow-extra-chr
plink --bfile outputs/genes --het --out outputs/genes --allow-extra-chr
plink --bfile outputs/nogea --het --out outputs/nogea --allow-extra-chr
plink --bfile outputs/nogea_thinned --het --out outputs/nogea_thinned --allow-extra-chr

# Calculate pairwise genetic distance
plink --bfile outputs/gea --allow-extra-chr --autosome-num 95 --distance square 1-ibs --const-fid --out outputs/gea 
plink --bfile outputs/genes --allow-extra-chr --autosome-num 95 --distance square 1-ibs --const-fid --out outputs/genes
plink --bfile outputs/nogea --allow-extra-chr --autosome-num 95 --distance square 1-ibs --const-fid --out outputs/nogea

# Run PCA
plink --bfile outputs/gea --allow-extra-chr --autosome-num 95 --pca 3 --out outputs/gea
plink --bfile outputs/genes --allow-extra-chr --autosome-num 95 --pca 3 --out outputs/genes
plink --bfile outputs/nogea --allow-extra-chr --autosome-num 95 --pca 3 --out outputs/nogea

# WINDOWED PI -------------------------------------------------------------------------
VCF=../../data/ccgp_data/58-Sceloporus_complete_coords_annotated.vcf.gz
vcftools --gzvcf $VCF --window-pi 10000 --out outputs/58-Sceloporus_10kb_windowpi
vcftools --gzvcf $VCF --window-pi 100000 --out outputs/58-Sceloporus_100kb_windowpi

# WINGEN FILES -------------------------------------------------------------------------
# Create dosage files
plink --bfile outputs/genes --recode A --allow-extra-chr --out outputs/genes # creates genes.raw
plink --bfile outputs/nogea_thinned --recode A --allow-extra-chr --out outputs/nogea_thinned # creates nogea_thinned.raw