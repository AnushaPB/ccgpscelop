conda activate ccgpscelop

# Genome and VCF files
GENOME=../../data/genome/rSceOcc1.20250428.p_ctg_chr_names.fna

# Using the VCF with all scaffolds
VCF=../../data/ccgp_data/58-Sceloporus_complete_coords_annotated_occidentalis_only.vcf.gz

# Annotation GFF converted using p1_gff_conversion.R
ANNOTATION=outputs/converted_annotation.gff

# Use csq to classify mutations
# Since this is unphased data, we use --local-csq to analyze each variant independently
bcftools csq -f "$GENOME" -g "$ANNOTATION" --local-csq "$VCF" -Ob -o outputs/csq.bcf
wc -l outputs/csq.bcf

# Extract all fields (no filtering) to a file
bcftools query -f '%CHROM\t%POS\t%INFO/BCSQ\n' outputs/csq.bcf > outputs/all_variant_csq.txt
wc -l outputs/all_variant_csq.txt