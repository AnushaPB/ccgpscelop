# run in analysis/gea
# Create bed file that has the structure contig snp_position snp_position (tab separated) from the GEA results
# Run in analysis/gea
# Intersect with the annotated .gff
conda activate ccgpscelop
GFF=../../data/annotated_genome/annotation/protein-coding.w_func_relabelled.gff
GEA=outputs/rda_sig_p01.bed
bedtools intersect -wo -a $GEA -b $GFF > outputs/gea_annotated.bed

# Grep out the GEA genes
grep -P "\tgene\t" outputs/gea_annotated.bed > outputs/gea_genes.bed

# THINGS TO CHECK
# - how many exons in each gene (don't trust if <2)
# - including linked SNPs in list