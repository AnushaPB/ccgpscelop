# Run in analysis/gea
source activate ccgpscelop

# Get plink 
PLINK=../../data/ccgp_data/58-Sceloporus_complete_coords_annotated
RDASNPS=outputs/bio1ndvi_rda_ids.txt
GENESNPS=outputs/bio1ndvi_gene_ids.txt

# Subset plink file with RDA snps (created using p1_process_rda.R)
plink --bfile $PLINK \
      --extract $RDASNPS \
      --make-bed \
      --out outputs/gea --allow-extra-chr

# Subset plink file with RDA genes (created using p3_intersect_genes.R)
plink --bfile outputs/gea \
      --extract $GENESNPS \
      --make-bed \
      --out outputs/genes --allow-extra-chr

# Subset plink file without RDA snps
plink --bfile $PLINK \
      --exclude $RDASNPS \
      --make-bed \
      --out outputs/nogea --allow-extra-chr

# Calculate allele frequencies for gea
plink --bfile outputs/gea --freq --out outputs/gea_freq --allow-extra-chr