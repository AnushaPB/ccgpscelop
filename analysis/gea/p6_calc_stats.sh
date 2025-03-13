# Run in analysis/gea
source activate ccgpscelop

# Get plink 
PLINK=../../data/ccgp_data/58-Sceloporus_complete_coords_annotated
VCF=../../data/ccgp_data/58-Sceloporus_complete_coords_annotated.vcf.gz
RDASNPS=outputs/bio1ndvi_rda_ids.txt
GENESNPS=outputs/bio1ndvi_gea_gene_ids.txt
SYNGENESNPS=outputs/bio1ndvi_gea_gene_syn_ids.txt
NONSYNGENESNPS=outputs/bio1ndvi_gea_gene_nonsyn_ids.txt

# Subset plink file with RDA snps (created using p1_process_rda.R)
plink --bfile $PLINK \
      --extract $RDASNPS \
      --make-bed \
      --out outputs/gea --allow-extra-chr

# Subset plink file with RDA genes (created using p5_intersect_genes.R)
plink --bfile outputs/gea \
      --extract $GENESNPS \
      --make-bed \
      --out outputs/genes --allow-extra-chr

# Subset plink file without RDA snps
plink --bfile $PLINK \
      --exclude $RDASNPS \
      --make-bed \
      --out outputs/nogea --allow-extra-chr

# Subset plink file to include RDA snps but not SNPs in genes
plink --bfile outputs/gea \
      --exclude $GENESNPS \
      --make-bed \
      --out outputs/geanogenes --allow-extra-chr

# Subset nonsynonymous GEA variants in genes
plink --bfile outputs/gea \
      --extract $NONSYNGENESNPS \
      --make-bed \
      --out outputs/nonsyn --allow-extra-chr

# Subset synonymous GEA varients in genes
plink --bfile outputs/gea \
      --extract $SYNGENESNPS \
      --make-bed \
      --out outputs/syn --allow-extra-chr

# Subset nogea plink file to the same number of snps as genes plink file
plink --bfile outputs/genes --write-snplist --out outputs/genes_snplist --allow-extra-chr
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
plink --bfile outputs/geanogenes --het --out outputs/geanogenes --allow-extra-chr
plink --bfile outputs/nogea_thinned --het --out outputs/nogea_thinned --allow-extra-chr
plink --bfile outputs/nonsyn --het --out outputs/nonsyn --allow-extra-chr
plink --bfile outputs/syn --het --out outputs/syn --allow-extra-chr

# WINDOWED PI -------------------------------------------------------------------------
vcftools --gzvcf $VCF --window-pi 10000 --out outputs/58-Sceloporus_10kb_windowpi
vcftools --gzvcf $VCF --window-pi 100000 --out outputs/58-Sceloporus_100kb_windowpi

# FST AND TAJIMA D  --------------------------------------------------

# FOR S. CA (POP 9)
vcftools --gzvcf $VCF --keep ../admixture/outputs/k9_pop9.txt --TajimaD 10000 --out outputs/58-Sceloporus_10kb_tajimad_pop9
vcftools --gzvcf $VCF --keep ../admixture/outputs/k9_pop9.txt --TajimaD 50000 --out outputs/58-Sceloporus_50kb_tajimad_pop9

# Coarse scale Tajima's D
for POP in {1..9}; do
  vcftools --gzvcf $VCF --keep ../admixture/outputs/k9_pop${POP}.txt --TajimaD 50000 --out outputs/58-Sceloporus_50kb_tajimad_pop${POP}
done

# Fine scale Tajima's D
for POP in {1..9}; do
  vcftools --gzvcf $VCF --keep ../admixture/outputs/k9_pop${POP}.txt --TajimaD 10000 --out outputs/58-Sceloporus_10kb_tajimad_pop${POP}
done

vcftools --gzvcf $VCF --weir-fst-pop ../admixture/outputs/k9_pop9.txt --weir-fst-pop ../admixture/outputs/k9_pop6.txt --fst-window-size 50000 --out outputs/58-Sceloporus_50kb_fst_pop9pop6

# WINGEN FILES -------------------------------------------------------------------------
# Create dosage files
plink --bfile outputs/genes --recode A --allow-extra-chr --out outputs/genes # creates genes.raw
plink --bfile outputs/nogea_thinned --recode A --allow-extra-chr --out outputs/nogea_thinned # creates nogea_thinned.raw


# PCA ----------------------------------------------------------------------------------
plink --bfile outputs/gea --allow-extra-chr --autosome-num 95 --pca 3 --out outputs/gea
plink --bfile outputs/genes --allow-extra-chr --autosome-num 95 --pca 3 --out outputs/genes
plink --bfile outputs/nogea --allow-extra-chr --autosome-num 95 --pca 3 --out outputs/nogea

# NOT SURE IF USING?

# Calculate pairwise genetic distance
plink --bfile outputs/gea --allow-extra-chr --autosome-num 95 --distance square 1-ibs --const-fid --out outputs/gea 
plink --bfile outputs/genes --allow-extra-chr --autosome-num 95 --distance square 1-ibs --const-fid --out outputs/genes
plink --bfile outputs/nogea --allow-extra-chr --autosome-num 95 --distance square 1-ibs --const-fid --out outputs/nogea
