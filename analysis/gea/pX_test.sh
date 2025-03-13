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

# Get min frequency of RDA SNPs
# gea_freq.frq
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.003472 0.099380 0.213300 0.231626 0.356200 0.500000 

# Filter nogea to only include SNPs with a min frequency of 0.1 and a max frequency of 0.4

plink --bfile outputs/nogea \
      --maf 0.1 \
      --max-maf 0.4 \
      --make-bed \
      --out outputs/nogea_frqfilter --allow-extra-chr

# Calculate heterozygosity
plink --bfile outputs/nogea_frqfilter --het --out outputs/nogea_frqfilter --allow-extra-chr


plink --bfile $PLINK \
      --thin-count 500000 \
      --make-bed \
      --out outputs/thinned --allow-extra-chr
plink --bfile outputs/thinned \
      --recode A \
      --out outputs/thinned --allow-extra-chr