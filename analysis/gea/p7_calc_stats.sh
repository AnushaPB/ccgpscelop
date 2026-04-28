# Run in analysis/gea
source activate ccgpscelop
# Get plink 
PLINK=../../data/ccgp_data/58-Sceloporus_complete_coords_annotated
VCF=../../data/ccgp_data/58-Sceloporus_complete_coords_annotated.vcf.gz
RDASNPS=outputs/bio1ndvi_rda_ids.txt
GENESNPS=outputs/bio1ndvi_gea_genes_ids.txt
NONSYNGENESNPS=outputs/bio1ndvi_gea_genes_nonsyn_ids.txt
ALLGENES=outputs/all_genes.bed

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

# Subset nonsynonymous GEA variants in genes
plink --bfile outputs/gea \
      --extract $NONSYNGENESNPS \
      --make-bed \
      --out outputs/nonsyn --allow-extra-chr

# Calculate allele frequencies
plink --bfile outputs/nonsyn \
      --freq --allow-extra-chr \
      --out outputs/nonsyn

# MAPPING FILES -------------------------------------------------------------------------
# Create dosage files for wingen and for gene maps
plink --bfile outputs/nonsyn --recode A --allow-extra-chr --out outputs/nonsyn # creates nonsyn.raw

# GENETIC DISTANCE MATRICES -------------------------------------------------------------
plink --bfile outputs/nonsyn --allow-extra-chr --autosome-num 95 --distance square 1-ibs --const-fid --out outputs/nonsyn