conda activate ccgpscelop
VCF=../../data/ccgp_data/58-Sceloporus_complete_coords_annotated.vcf.gz
SAMPLEIDS=../../data/final_sampleids.txt

vcftools --gzvcf $VCF --keep $SAMPLEIDS --window-pi 10000 --out outputs/58-Sceloporus_10kb_windowpi