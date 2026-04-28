conda activate ccgpscelop
VCF=../../data/ccgp_data/58-Sceloporus_complete_coords_annotated_final_samples.vcf.gz

# Index VCF file if not already indexed
if [ ! -f "${VCF}.tbi" ]; then
    tabix -p vcf $VCF
fi

# Non-overlapping windows of 10kb and 50kb
vcftools --gzvcf $VCF --window-pi 10000 --out outputs/58-Sceloporus_10kb_windowpi
vcftools --gzvcf $VCF --window-pi 50000 --out outputs/58-Sceloporus_50kb_windowpi

# Calculate within hot and cold regions
HOTIDS=../genetic_diversity/outputs/hot_individuals.txt
COLDIDS=../genetic_diversity/outputs/cold_individuals.txt
vcftools --gzvcf $VCF --keep $HOTIDS --window-pi 50000 --out outputs/hot_50kb_windowpi
vcftools --gzvcf $VCF --keep $COLDIDS --window-pi 50000 --out outputs/cold_50kb_windowpi

# Calculate pi within pops
POPDIR=../admixture/outputs
for POPFILE in ${POPDIR}/k8_pop{1..8}.txt; do
    POP=$(basename "$POPFILE" .txt)
    vcftools --gzvcf "$VCF" \
        --keep "$POPFILE" \
        --window-pi 50000 \
        --out "outputs/${POP}_50kb_windowpi"
done