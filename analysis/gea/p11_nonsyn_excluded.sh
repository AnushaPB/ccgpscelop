
conda activate ccgpscelop
vcf="../../data/ccgp_data/58-Sceloporus_complete_coords_annotated_final_samples.vcf.gz"
bed="outputs/all_nonsynonymous.bed"
out="outputs/58-Sceloporus_no_nonsyn_50kb"

# contig lengths from VCF header
bcftools view -h "$vcf" | \
awk -F'[=,>]' '/^##contig=<ID=/{print $3"\t"$5}' > outputs/genome.txt

# make 50 kb exclusion regions around nonsynonymous sites
bedtools slop -i "$bed" -g outputs/genome.txt -b 50000 | \
sort -k1,1 -k2,2n | \
bedtools merge -i - > outputs/nonsyn_50kb_merged.bed

# remove masked regions
bcftools view \
-T ^outputs/nonsyn_50kb_merged.bed \
-Oz \
-o ${out}.vcf.gz \
"$vcf"

bcftools index -t ${out}.vcf.gz

# calculate heterozygosity
plink2 \
--vcf ${out}.vcf.gz \
--double-id \
--allow-extra-chr \
--het \
--out ${out}