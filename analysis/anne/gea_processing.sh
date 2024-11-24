# Run in analysis/anne

# RUN PCA ON FULL DATASET

plink --vcf 58-Sceloporus_annotated_pruned_0.6.vcf.gz --pca 10 --out 58-Sceloporus_annotated_pruned_0.6 --allow-extra-chr --autosome-num 95 --const-fid

# EXTRACT ONE SCAFFOLD -------------------------------------------------------------------

# Extract a single scaffold
# CONTIG="Scaffold_10__1_contigs__length_2370820"
CONTIG="Scaffold_161__1_contigs__length_1846707"
vcftools --gzvcf 58-Sceloporus_annotated_pruned_0.6.vcf.gz --chr $CONTIG --recode --recode-INFO-all --out "58-Sceloporus_annotated_pruned_0.6_onecontig"

# RUN GEA --------------------------------------------------------------------------------