# run in analysis/gea
# Create bed file that has the structure contig snp_position snp_position (tab separated) from the GEA results
# Run in analysis/geabed > outputs/gea_genes.bed

# THINGS TO CHECK
# - how many exons in each gene (don't trust if <2)
# - including linked SNPs in list
# Intersect with the annotated .gff
conda activate ccgpscelop
#conda install -c bioconda bedops

GFF=../../data/annotated_genome/annotation/protein-coding.w_func_relabelled.gff

# Grep out all genes from the annotation
grep -P "\tgene\t" $GFF > outputs/all_genes.gff

# Convert to bed
gff2bed < outputs/all_genes.gff > outputs/all_genes.bed

# Note: make sure start and end are 1 SNP apart in created bed file (bed is 0-based)
head -n 5 outputs/bio1ndvi_gea.bed 

# Get the GEA genes from the environmental PCs
# -wo = write out the original A and B entries
bedtools intersect -wo -a outputs/bio1ndvi_gea.bed -b outputs/all_genes.bed > outputs/bio1ndvi_gea_genes.bed
#bedtools intersect -a outputs/pca_gea.bed -b outputs/all_genes.bed > outputs/pca_gea_genes.bed

# Check number of SNPs in genes
wc -l outputs/bio1ndvi_gea_genes.bed #631871 
