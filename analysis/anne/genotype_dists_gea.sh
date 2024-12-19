### Plotting genotypes of topmost outliers across landscape
# Run within analysis/anne directory
mamba activate ccgpscelop

tabix -p vcf ../../data/ccgp_data/58-Sceloporus_annotated_pruned_0.6.vcf.gz
bcftools view ../../data/ccgp_data/58-Sceloporus_annotated_pruned_0.6.vcf.gz --regions-file outputs/supersig_SNPs.txt > outputs/58-Sceloporus_chr1_supersigSNPs.vcf