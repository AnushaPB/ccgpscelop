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
bcftools view -i 'INFO/BCSQ ~ "^missense\|" | INFO/BCSQ ~ "^frameshift\|" | INFO/BCSQ ~ "^stop_gained\|" | INFO/BCSQ ~ "^stop_lost\|" | INFO/BCSQ ~ "^start_lost\|" | INFO/BCSQ ~ "^inframe_altering\|" | INFO/BCSQ ~ "^inframe_deletion\|" | INFO/BCSQ ~ "^inframe_insertion\|"' outputs/csq.bcf -Ov -o outputs/nonsynonymous.vcf

bcftools query -f '%CHROM\t%POS\t%INFO/BCSQ\n' outputs/nonsynonymous.vcf | \
awk -F'\t' '{ split($3, a, "\\|"); print $1 ":" $2 "\t" a[1] }' \
> outputs/nonsynonymous_variants.txt


# Count the number of non-synonymous variants
bcftools view -v snps outputs/nonsynonymous.vcf | grep -v "^#" | wc -l

# Extract the relevant information from the VCF file
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/BCSQ\n' outputs/nonsynonymous.vcf | \
awk -F'\t' '$5 ~ /^(missense|frameshift|stop_gained|stop_lost|start_lost|inframe_altering|inframe_deletion|inframe_insertion)\|/' \
> outputs/nonsynonymous_variants.txt

# Calculate length of txt file
wc -l outputs/nonsynonymous_variants.txt

# Remove the temporary file
rm fixed_annotation.gff3