# Run in analysis/gea
source activate ccgpscelop

# Get plink 
PLINK=../../data/ccgp_data/58-Sceloporus_complete_coords_annotated

# Files created from analysis/gea/p7_csq_intersect.R
ALLNONSYN=../gea/outputs/all_nonsynonymous.bed
ALLSYN=../gea/outputs/all_synonymous.bed
NONSYNGENESNPS=../gea/outputs/bio1ndvi_gea_genes_nonsyn_ids.txt

# Subset plink file with all nonsynonymous variants
plink --bfile $PLINK \
      --extract range $ALLNONSYN \
      --make-bed \
      --out outputs/allnonsyn --allow-extra-chr

# Subset plink file with all synonymous variants
plink --bfile $PLINK \
      --extract range $ALLSYN \
      --make-bed \
      --out outputs/allsyn --allow-extra-chr

# Subset nonsynonymous GEA variants in genes
plink --bfile $PLINK \
      --extract $NONSYNGENESNPS \
      --make-bed \
      --out outputs/nonsyn --allow-extra-chr
  
# Calculate heterozygosity stats
plink --bfile outputs/allnonsyn --het --out outputs/allnonsyn --allow-extra-chr
plink --bfile outputs/allsyn --het --out outputs/allsyn --allow-extra-chr
plink --bfile outputs/nonsyn --het --out outputs/nonsyn --allow-extra-chr

