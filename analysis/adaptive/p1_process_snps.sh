# Run in analysis/adaptive
mamba activate ccgpscelop

# ================= BIO1 + NDVI results

# Extract relevant SNPs from vcf using the bed file
bcftools view -R ../gea/outputs/bio1ndvi_gea_genes_nonsyn.bed ../../data/ccgp_data/58-Sceloporus_complete_coords_annotated.vcf.gz > ./outputs/58-Sceloporus_bio1ndvi_gea_genes_nonsyn.vcf
bcftools query -f "%POS\n" ./outputs/58-Sceloporus_bio1ndvi_gea_genes_nonsyn.vcf | wc -l # 29,445
cat ../gea/outputs/bio1ndvi_gea_genes_nonsyn.bed | wc -l # 29,426

bcftools view -R ../gea/outputs/bio1ndvi_gea.bed ../../data/ccgp_data/58-Sceloporus_complete_coords_annotated.vcf.gz > ./outputs/58-Sceloporus_bio1ndvi_gea.vcf
bcftools query -f "%POS\n" ./outputs/58-Sceloporus_bio1ndvi_gea.vcf | wc -l # 1,544,696
cat ../gea/outputs/bio1ndvi_gea.bed | wc -l # 1,543,125

# Now calculate Plink distances for the new vcf for running GDM
# plink --vcf ./outputs/58-Sceloporus_bio1ndvi_gea.vcf --out ./outputs/58-Sceloporus_bio1ndvi_gea --distance square --const-fid --allow-extra-chr
plink --vcf ./outputs/58-Sceloporus_bio1ndvi_gea.vcf --out ./outputs/58-Sceloporus_bio1ndvi_gea_ibs --distance square 1-ibs --const-fid --allow-extra-chr

# ================= rasterPCA results

# Take first two cols of gea_gene_snp.csv as the regions file
# awk '{print $1,$2 > "./outputs/pca_gea_gene_regions.csv"}' ../gea/outputs/pca_gea_genes.bed

# Convert to tsv and remove csv
# tr ',' '\t' < ./outputs/pca_gea_gene_regions.csv > ./outputs/pca_gea_gene_regions.txt
# rm ./outputs/pca_gea_gene_regions.csv

# Extract relevant SNPs from vcf
# bcftools view ../../data/ccgp_data/58-Sceloporus_complete_coords_annotated.vcf.gz --regions-file ./outputs/pca_gea_gene_regions.txt > ./outputs/58-Sceloporus_pca_gea_genes.vcf
# bcftools query -f "%POS\n" ./outputs/58-Sceloporus_pca_gea_genes.vcf | wc -l