# Run in analysis/adaptive
mamba activate ccgpscelop

# ================= Additional processing of vcf file

# Remove irrelevant individuals from vcf (resulting n=257 inds)
bcftools view -s "^Scelocci_CAS213197,Scelocci_CAS214858" ../../data/ccgp_data/58-Sceloporus_complete_coords_annotated_occidentalis_only.vcf.gz > ./outputs/58-Sceloporus_complete_coords_annotated_occidentalis_only_correct.vcf
# Retain autosomes only
bgzip ./outputs/58-Sceloporus_complete_coords_annotated_occidentalis_only_correct.vcf
tabix -p vcf ./outputs/58-Sceloporus_complete_coords_annotated_occidentalis_only_correct.vcf.gz
bcftools index --stats ./outputs/58-Sceloporus_complete_coords_annotated_occidentalis_only_correct.vcf.gz | awk '{print $1}' > ./outputs/chrom_list # list of regions
# Remove two sex chromosomes and save as autosome_list

# Format properly for bcftools & only retain only autosomes
tr '\n' ',' < ./outputs/autosome_list > ./outputs/autosome_list_formatted # replace newlines with commas
bcftools view ./outputs/58-Sceloporus_complete_coords_annotated_occidentalis_only_correct.vcf.gz --regions `cat ./outputs/autosome_list_formatted` > ./outputs/58-Sceloporus_complete_coords_annotated_occidentalis_only_correct_autosomes.vcf
# Check that there are 35 chroms in final vcf
bcftools query -f '%CHROM\n' ./outputs/58-Sceloporus_complete_coords_annotated_occidentalis_only_correct_autosomes.vcf | sort | uniq | wc -l
bgzip ./outputs/58-Sceloporus_complete_coords_annotated_occidentalis_only_correct_autosomes.vcf
tabix -p vcf ./outputs/58-Sceloporus_complete_coords_annotated_occidentalis_only_correct_autosomes.vcf.gz

# ================= BIO1 + NDVI results

# Extract relevant SNPs from vcf using the bed file
bcftools view -R ../gea/outputs/bio1ndvi_gea_genes_nonsyn.bed ./outputs/58-Sceloporus_complete_coords_annotated_occidentalis_only_correct_autosomes.vcf.gz > ./outputs/58-Sceloporus_bio1ndvi_gea_genes_nonsyn.vcf

# Now calculate Plink distances for the new vcf for running GDM
plink --vcf ./outputs/58-Sceloporus_bio1ndvi_gea_genes_nonsyn.vcf --out ./outputs/58-Sceloporus_bio1ndvi_gea_ibs --distance square 1-ibs --const-fid --allow-extra-chr

# ================= rasterPCA results

# Take first two cols of gea_gene_snp.csv as the regions file
# awk '{print $1,$2 > "./outputs/pca_gea_gene_regions.csv"}' ../gea/outputs/pca_gea_genes.bed

# Convert to tsv and remove csv
# tr ',' '\t' < ./outputs/pca_gea_gene_regions.csv > ./outputs/pca_gea_gene_regions.txt
# rm ./outputs/pca_gea_gene_regions.csv

# Extract relevant SNPs from vcf
# bcftools view ../../data/ccgp_data/58-Sceloporus_complete_coords_annotated.vcf.gz --regions-file ./outputs/pca_gea_gene_regions.txt > ./outputs/58-Sceloporus_pca_gea_genes.vcf
# bcftools query -f "%POS\n" ./outputs/58-Sceloporus_pca_gea_genes.vcf | wc -l

# ================= STATS TO REPORT

# Number of chroms in original data
cat ./outputs/chrom_list | wc -l # 37
# Number of chroms in analyzed data (i.e., autosomes only)
cat ./outputs/autosome_list | wc -l # 35
# Number of sites in bed file
cat ../gea/outputs/bio1ndvi_gea_genes_nonsyn.bed | wc -l # 8860 variants
# Verify that the number of sites in analyzed vcf match the bed file
bcftools query -f "%POS\n" ./outputs/58-Sceloporus_bio1ndvi_gea_genes_nonsyn.vcf | wc -l # 8860
# Number of individuals in analyzed vcf
bcftools query -l ./outputs/58-Sceloporus_bio1ndvi_gea_genes_nonsyn.vcf | wc -l # 257