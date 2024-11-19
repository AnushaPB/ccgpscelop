
# Get plink 
PLINK=../../data/ccgp_data/58-Sceloporus_annotated_pruned_0.6

# Subset plink file with RDA snps (created using process_rda.R)
plink --bfile $PLINK \
      --extract outputs/rda_ids.txt \
      --make-bed \
      --out outputs/gea --allow-extra-chr

# Subset plink file with genes (created using intersect_genes.R)
plink --bfile outputs/gea \
      --extract outputs/gene_ids.txt \
      --make-bed \
      --out outputs/genes --allow-extra-chr

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

# Calculate pairwise genetic distance
plink --bfile outputs/gea --allow-extra-chr --autosome-num 95 --distance square 1-ibs --const-fid --out outputs/gea 
plink --bfile outputs/genes --allow-extra-chr --autosome-num 95 --distance square 1-ibs --const-fid --out outputs/genes 
plink --bfile $PLINK --allow-extra-chr --autosome-num 95 --distance square 1-ibs --const-fid --out outputs/neutral 

# Run PCA
plink --bfile outputs/gea --allow-extra-chr --autosome-num 95 --pca 3 --out outputs/gea
plink --bfile outputs/genes --allow-extra-chr --autosome-num 95 --pca 3 --out outputs/genes
plink --bfile $PLINK --allow-extra-chr --autosome-num 95 --pca 3 --out outputs/neutral