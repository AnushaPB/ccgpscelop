GFF=../../data/spea_contaminated_genome/contaminated_annotated_genome/annotated_genome/annotation/protein-coding.w_func_relabelled.gff


grep -P "\tgene\t" $GFF > outputs/all_genes_old.gff

wc -l outputs/all_genes_old.gff #74,808 genes

# Compariosn with new genes
grep -P "\tgene\t" $GFF > outputs/all_genes.gff
gff2bed < outputs/all_genes.gff > outputs/all_genes.bed

wc -l outputs/all_genes.gff #20,638 genes

gff2bed < outputs/all_genes_old.gff > outputs/all_genes_old.bed