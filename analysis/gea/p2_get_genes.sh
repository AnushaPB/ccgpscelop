# run in analysis/gea
# Create bed file that has the structure contig snp_position snp_position (tab separated) from the GEA results
# Run in analysis/gea
# Intersect with the annotated .gff
conda activate ccgpscelop
GFF=../../data/annotated_genome/annotation/protein-coding.w_func_relabelled.gff
PCA=outputs/rda_sig_p01.bed
BIO1=outputs/bio1_sig.bed

# Grep out all genes from the annotation
grep -P "\tgene\t" $GFF > outputs/all_genes.bed

# Get the GEA genes from the environmental PCs
bedtools intersect -wo -a $PCA -b outputs/all_genes.bed > outputs/gea_genes.bed

# Get the GEA genes from BIO1 
bedtools intersect -wo -a $BIO1 -b outputs/all_genes.bed > outputs/bio1_genes.bed

# THINGS TO CHECK
# - how many exons in each gene (don't trust if <2)
# - including linked SNPs in list