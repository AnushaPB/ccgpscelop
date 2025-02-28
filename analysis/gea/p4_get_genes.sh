# run in analysis/gea
# Create bed file that has the structure contig snp_position snp_position (tab separated) from the GEA results
# Run in analysis/geabed > outputs/gea_genes.bed

# THINGS TO CHECK
# - how many exons in each gene (don't trust if <2)
# - including linked SNPs in list
# Intersect with the annotated .gff
conda activate ccgpscelop
GFF=../../data/annotated_genome/annotation/protein-coding.w_func_relabelled.gff

# Grep out all genes from the annotation
grep -P "\tgene\t" $GFF > outputs/all_genes.bed

# Get the GEA genes from the environmental PCs
bedtools intersect -wo -a outputs/bio1ndvi_gea.bed -b outputs/all_genes.bed > outputs/bio1ndvi_gea_genes.bed
#bedtools intersect -wo -a outputs/pca_gea.bed -b outputs/all_genes.bed > outputs/pca_gea_genes.bed