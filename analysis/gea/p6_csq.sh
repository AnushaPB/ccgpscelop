conda activate ccgpscelop

GENOME=../../data/annotated_genome/jordan-uni4378-mb-hirise-65qvq__12-23-2023__final_assembly_relabelled.fasta
ANNOTATION=../../data/annotated_genome/annotation/protein-coding.w_func_relabelled.gff
VCF=../../data/ccgp_data/58-Sceloporus_complete_coords_annotated.vcf.gz

# Convert GFF to csq format (CHECK THIS)
python convert_gff.py $ANNOTATION outputs/fixed_annotation.gff

# Use csq to classify mutations
bcftools csq -f "$GENOME" -g outputs/fixed_annotation.gff --local-csq "$VCF" -Ob -o outputs/csq.bcf

# Extract all fields (no filtering) to a file
bcftools query -f '%CHROM\t%POS\t%INFO/BCSQ\n' outputs/csq.bcf > outputs/all_variant_csq.txt
wc -l outputs/all_variant_csq.txt

# Remove the temporary file
rm outputs/fixed_annotation.gff