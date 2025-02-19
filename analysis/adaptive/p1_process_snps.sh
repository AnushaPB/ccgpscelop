# ================= rasterPCA results
# Take first two cols of gea_gene_snp.csv as the regions file
# awk '{print $1,$2 > "./outputs/pca_gea_gene_regions.csv"}' ../gea/outputs/pca_gea_genes.bed

# Convert to tsv and remove csv
# tr ',' '\t' < ./outputs/pca_gea_gene_regions.csv > ./outputs/pca_gea_gene_regions.txt
# rm ./outputs/pca_gea_gene_regions.csv

# Extract relevant SNPs from vcf
# bcftools view ../../data/ccgp_data/58-Sceloporus_complete_coords_annotated.vcf.gz --regions-file ./outputs/pca_gea_gene_regions.txt > ./outputs/58-Sceloporus_pca_gea_genes.vcf
# bcftools query -f "%POS\n" ./outputs/58-Sceloporus_pca_gea_genes.vcf | wc -l

# ================= BIO1 + NDVI results
# BIO1 + NDVI outliers
awk '{print $1,$2 > "./outputs/bio1ndvi_gea_gene_regions.csv"}' ../gea/outputs/bio1ndvi_gea_genes.bed

# Convert to tsv and remove csv
tr ' ' '\t' < ./outputs/bio1ndvi_gea_gene_regions.csv > ./outputs/bio1ndvi_gea_gene_regions.txt
rm ./outputs/bio1ndvi_gea_gene_regions.csv

# Extract relevant SNPs from vcf
bcftools view ../../data/ccgp_data/58-Sceloporus_complete_coords_annotated.vcf.gz --regions-file ./outputs/bio1ndvi_gea_gene_regions.txt > ./outputs/58-Sceloporus_bio1ndvi_gea_genes.vcf
# bcftools query -f "%POS\n" ./outputs/58-Sceloporus_bio1ndvi_gea_genes.vcf | wc -l # 600056 vs 63294