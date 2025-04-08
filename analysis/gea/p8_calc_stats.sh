# Run in analysis/gea
source activate ccgpscelop

# Get plink 
PLINK=../../data/ccgp_data/58-Sceloporus_complete_coords_annotated
VCF=../../data/ccgp_data/58-Sceloporus_complete_coords_annotated.vcf.gz
RDASNPS=outputs/bio1ndvi_rda_ids.txt
GENESNPS=outputs/bio1ndvi_gea_gene_ids.txt
#SYNGENESNPS=outputs/bio1ndvi_gea_gene_syn_ids.txt
SYNGENESNPS=outputs/bio1ndvi_gea_gene_unlinked_syn_ids.txt
NONSYNGENESNPS=outputs/bio1ndvi_gea_gene_nonsyn_ids.txt
NONSYNNOTGEAGENESNPS=outputs/notgeagenes_nonsyn.bed
ALLGENES=outputs/all_genes.bed
ALLNONSYN=outputs/all_nonsynonymous.bed

# Add buffer
awk -v buf=10000 'BEGIN{OFS="\t"} {start=$2-buf; end=$3+buf; if(start<0) start=0; print $1, start, end}' $ALLGENES > outputs/all_genes_buffer10kb.bed
# Format as range file
awk '{print $1, $2, $3, "region"NR}' OFS="\t" outputs/all_genes_buffer10kb > outputs/all_genes_buffer10kb.bed
ALLGENESBUFFERED=outputs/all_genes_buffer10kb.bed

# Subset plink file with GEA snps (created using p1_process_rda.R)
plink --bfile $PLINK \
      --extract $RDASNPS \
      --make-bed \
      --out outputs/gea --allow-extra-chr

# Subset plink file with GEA genes (created using p5_intersect_genes.R)
plink --bfile outputs/gea \
      --extract $GENESNPS \
      --make-bed \
      --out outputs/genes --allow-extra-chr

# Subset plink file without GEA snps
plink --bfile $PLINK \
      --exclude $RDASNPS \
      --make-bed \
      --out outputs/nogea --allow-extra-chr

# Subset plink file with no GEA SNPs and no genes + buffer
plink --bfile outputs/nogea \
      --exclude range $ALLGENESBUFFERED \
      --make-bed \
      --out outputs/nogeanogenes --allow-extra-chr

# Subset nonsynonymous GEA variants in genes
plink --bfile outputs/gea \
      --extract $NONSYNGENESNPS \
      --make-bed \
      --out outputs/nonsyn --allow-extra-chr

# Subset synonymous GEA variants in genes
plink --bfile outputs/gea \
      --extract $SYNGENESNPS \
      --make-bed \
      --out outputs/syn --allow-extra-chr

# Subset non-synonymous variants in genes not found by GEA
# NEED TO ADD REGION COLUMN FOR PLINK (DOESN'T MATTER WHAT IT IS)
awk '{print $1, $2, $3, "region"NR}' OFS="\t" $NONSYNNOTGEAGENESNPS > outputs/notgeagenes_nonsyn_fixed.bed
plink --bfile $PLINK \
      --extract range  outputs/notgeagenes_nonsyn_fixed.bed \
      --make-bed \
      --out outputs/nonsynnotgea --allow-extra-chr
rm outputs/notgeagenes_nonsyn_fixed.bed

# Get all non-synonymous SNPs
plink --bfile $PLINK \
      --extract range $ALLNONSYN \
      --make-bed \
      --out outputs/allnonsyn --allow-extra-chr
   
# Subset nogeanogenes plink file to the same number of snps as non-synonymous genes plink file
plink --bfile outputs/nonsyn --write-snplist --out outputs/nonsyn_snplist --allow-extra-chr
N=$(wc -l < outputs/nonsyn_snplist.snplist)
plink --bfile outputs/nogeanogenes --thin-count $N --out outputs/nogeanogenes_thinned --make-bed --allow-extra-chr

# Calculate allele frequencies
plink --bfile outputs/nogeanogenes \
      --freq --allow-extra-chr \
      --out outputs/nogeanogenes

plink --bfile outputs/nonsyn \
      --freq --allow-extra-chr \
      --out outputs/nonsyn

plink --bfile outputs/allnonsyn \
      --freq --allow-extra-chr \
      --out outputs/allnonsyn      

Rscript script_pick_frqs.R

plink --bfile outputs/nogeanogenes \
      --extract outputs/neutral_snp_frqmatch.txt \
      --make-bed \
      --out outputs/frqmatch --allow-extra-chr   

plink --bfile $PLINK \
      --extract outputs/allnonsyn_snp_frqmatch.txt \
      --make-bed \
      --out outputs/allnonsyn_frqmatch --allow-extra-chr   

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
plink --bfile outputs/nogeanogenes --het --out outputs/nogeanogenes --allow-extra-chr
plink --bfile outputs/nogeanogenes_thinned --het --out outputs/nogeanogenes_thinned --allow-extra-chr
plink --bfile outputs/nonsyn --het --out outputs/nonsyn --allow-extra-chr
plink --bfile outputs/syn --het --out outputs/syn --allow-extra-chr
plink --bfile outputs/frqmatch --het --out outputs/frqmatch --allow-extra-chr
plink --bfile outputs/allnonsyn_frqmatch --het --out outputs/allnonsyn_frqmatch --allow-extra-chr
plink --bfile outputs/nonsynnotgea --het --out outputs/nonsynnotgea --allow-extra-chr
plink --bfile outputs/allnonsyn --het --out outputs/allnonsyn --allow-extra-chr

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

# MAPPING FILES -------------------------------------------------------------------------
# Create dosage files for wingen and for gene maps
plink --bfile outputs/nonsyn --recode A --allow-extra-chr --out outputs/genes # creates genes.raw
plink --bfile outputs/nogea_thinned --recode A --allow-extra-chr --out outputs/nogea_thinned # creates nogea_thinned.raw

# PCA ----------------------------------------------------------------------------------
plink --bfile outputs/gea --allow-extra-chr --autosome-num 95 --pca 3 --out outputs/gea
plink --bfile outputs/nonsyn --allow-extra-chr --autosome-num 95 --pca 3 --out outputs/genes
plink --bfile outputs/nogea --allow-extra-chr --autosome-num 95 --pca 3 --out outputs/nogea

# NOT SURE IF USING?

# Calculate pairwise genetic distance
plink --bfile outputs/gea --allow-extra-chr --autosome-num 95 --distance square 1-ibs --const-fid --out outputs/gea 
plink --bfile outputs/nonsyn --allow-extra-chr --autosome-num 95 --distance square 1-ibs --const-fid --out outputs/genes
plink --bfile outputs/nogea --allow-extra-chr --autosome-num 95 --distance square 1-ibs --const-fid --out outputs/nogea
