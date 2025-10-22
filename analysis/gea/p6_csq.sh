conda activate ccgpscelop

# Install agat to handle GFF files if not already installed
#conda install -c bioconda agat

GENOME=../../data/genome/rSceOcc1.20250428.p_ctg_chr_names.fna
ANNOTATION=../../data/genome/annotation/complete.genomic.gff
VCF=../../data/ccgp_data/58-Sceloporus_complete_coords_annotated_occidentalis_only.vcf.gz

# Convert GFF to csq format (CHECK THIS)
python script_convert_gff.py $ANNOTATION outputs/fixed_annotation.gff

# wget https://github.com/samtools/bcftools/blob/e4a21b8273eb8c79a03ce6e3e651d5b5422ff75b/misc/gff2gff
# chmod +x gff2gff
# cat "$ANNOTATION" | ./gff2gff | bgzip -c > outputs/fixed_annotation.gff.gz
# zcat outputs/fixed_annotation.gff.gz | ./gff2gff | gzip -c > outputs/fixed_annotation2.gff.gz

# Use csq to classify mutations
# Since this is unphased data, we use --local-csq to analyze each variant independently
bcftools csq -f "$GENOME" -g outputs/fixed_annotation.gff --local-csq "$VCF" -Ob -o outputs/csq.bcf

# Extract all fields (no filtering) to a file
bcftools query -f '%CHROM\t%POS\t%INFO/BCSQ\n' outputs/csq.bcf > outputs/all_variant_csq.txt
wc -l outputs/all_variant_csq.txt

# Remove the temporary file
rm outputs/fixed_annotation.gff