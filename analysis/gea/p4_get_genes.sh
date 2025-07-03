# run in analysis/gea
# Create bed file that has the structure contig snp_position snp_position (tab separated) from the GEA results
# Run in analysis/geabed > outputs/gea_genes.bed

# THINGS TO CHECK
# - how many exons in each gene (don't trust if <2)
# - including linked SNPs in list
# Intersect with the annotated .gff
conda activate ccgpscelop
#conda install -c bioconda bedops

GFF=../../data/genome/annotation/complete.genomic.gff

# 1. Extract all gene entries and convert to BED
grep -P "\tgene\t" $GFF > outputs/all_genes.gff
gff2bed < outputs/all_genes.gff > outputs/all_genes.bed

# 2. Extract all exon (coding region + UTR) entries and convert to BED
grep -P "\texon\t" $GFF > outputs/all_exons.gff
gff2bed < outputs/all_exons.gff > outputs/all_exons.bed

# 3. Extract all CDS (coding region) entries and convert to BED
grep -P "\tCDS\t" $GFF > outputs/all_cds.gff
gff2bed < outputs/all_cds.gff > outputs/all_cds.bed

# Note: make sure start and end are 1 SNP apart in created bed file (bed is 0-based)
head -n 5 outputs/bio1ndvi_gea.bed 

# Get the GEA genes from the environmental PCs
# -wo = write out the original A and B entries
# "in genes" means SNPs fall within the coding region of the gene
bedtools intersect -wo -a outputs/bio1ndvi_gea.bed -b outputs/all_cds.bed > outputs/bio1ndvi_gea_genes.bed

# Check number of SNPs in genes
wc -l outputs/bio1ndvi_gea_genes.bed #34,410 
