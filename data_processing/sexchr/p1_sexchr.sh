# Window depth
VCF="../../data/ccgp_data/58-Sceloporus_complete_coords_annotated.vcf.gz"
OUTDIR="outputs"
mkdir -p "$OUTDIR"

# Genome-wide mean depth per sample (across *all* contigs)
vcftools --gzvcf "$VCF" \
  --depth \
  --out "$OUTDIR/genomewide"

sex_linked_scaffolds=("sex_linked_1" "sex_linked_2")
REGIONS=$(IFS=,; echo "${sex_linked_scaffolds[*]}")

# list samples in the same order as columns in the DP matrix
bcftools query -l "$VCF" > "$OUTDIR/samples.txt"

# Per-site DP table for the two scaffolds: CHROM POS then one DP column per sample
bcftools query -r "$REGIONS" \
  -f '%CHROM\t%POS[\t%DP]\n' "$VCF" \
  > "$OUTDIR/sex_linked.dp.tsv"