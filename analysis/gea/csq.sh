conda activate ccgpscelop

GENOME=../../data/annotated_genome/jordan-uni4378-mb-hirise-65qvq__12-23-2023__final_assembly_relabelled.fasta
ANNOTATION=../../data/annotated_genome/annotation/protein-coding.w_func_relabelled.gff
VCF=../../data/ccgp_data/58-Sceloporus_complete_coords_annotated.vcf.gz
 fixed_annotation.gff3

# Convert GFF to csq format (CHECK THIS)
python convert_gff.py $ANNOTATION fixed_annotation.gff

# Use csq to classify mutations
bcftools csq -f "$GENOME" -g fixed_annotation.gff --local-csq "$VCF" -Ob -o outputs/csq.bcf

# Filter all nonsynonymous variants (missense + other functional changes)
bcftools view -i 'INFO/BCSQ ~ "missense"' outputs/csq.bcf -Ov -o outputs/missense.vcf

bcftools view -i 'INFO/BCSQ ~ "missense|frameshift|inframe_altering|inframe_deletion|inframe_insertion|start_lost|stop_gained|stop_lost"' outputs/csq.bcf -Ov -o outputs/nonsynonymous.vcf

# Count the number of non-synonymous variants
bcftools view -v snps outputs/nonsynonymous.vcf | grep -v "^#" | wc -l
bcftools view -v snps outputs/missense.vcf | grep -v "^#" | wc -l

# Extract the relevant information from the VCF file
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/BCSQ\n' outputs/nonsynonymous.vcf | \
awk -F'\t' '$5 ~ /missense|frameshift|inframe_altering|inframe_deletion|inframe_insertion|start_lost|stop_gained|stop_lost/' > outputs/nonsynonymous_variants.txt

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/BCSQ\n' outputs/missense.vcf  | \
awk -F'\t' '$5 ~ /missense|frameshift|inframe_altering|inframe_deletion|inframe_insertion|start_lost|stop_gained|stop_lost/' > outputs/nonsynonymous_variants.txt

# Calculate length of txt file
wc -l outputs/nonsynonymous_variants.txt
# Remove the temporary file
rm fixed_annotation.gff3