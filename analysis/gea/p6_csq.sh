conda activate ccgpscelop

GENOME=../../data/genome/rSceOcc1.20250428.p_ctg_chr_names.fna
ANNOTATION=../../data/genome/annotation/complete.genomic.gff
VCF=../../data/ccgp_data/58-Sceloporus_complete_coords_annotated.vcf.gz

# Convert GFF to csq format (CHECK THIS)
python script_convert_gff.py $ANNOTATION outputs/fixed_annotation.gff

# Use csq to classify mutations
bcftools csq -f "$GENOME" -g outputs/fixed_annotation.gff --local-csq "$VCF" -Ob -o outputs/csq.bcf

# Extract all fields (no filtering) to a file
bcftools query -f '%CHROM\t%POS\t%INFO/BCSQ\n' outputs/csq.bcf > outputs/all_variant_csq.txt
wc -l outputs/all_variant_csq.txt

# Remove the temporary file
rm outputs/fixed_annotation.gff