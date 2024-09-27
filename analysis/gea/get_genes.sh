# Create bed file that has the structure contig snp_position snp_position (tab separated) from the GEA results

# Intersect with the annotated .gff
GFF=../../data/annotation/protein-coding.w_func_relabelled.gff
GEA=gea.bed
bedtoools intersect -wo -a $GEA -b $GFF > gea_annotated.bed

# Grep out the genes
grep -P "\tgene\t" gea_annotated.bed > gea_genes.bed

# THINGS TO CHECK
# - how many exons in each gene (don't trust if <2)
# - including linked SNPs in list