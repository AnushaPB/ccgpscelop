
# Calculate depth just for sex_linked_1 and sex_linked_2
VCF="../../data/ccgp_data/58-Sceloporus_complete_coords_annotated.vcf.gz"
PLINK="../../data/ccgp_data/58-Sceloporus_complete_coords_annotated"
sex_linked_scaffolds=("sex_linked_1" "sex_linked_2")

for scaffold in "${sex_linked_scaffolds[@]}"; do
  vcftools --gzvcf "$VCF" \
            --chr "$scaffold" \
            --depth \
            --out "outputs/${scaffold}"
done

# Compute heterozygosity 
plink --bfile $PLINK --het --out outputs/sex_linked_1 --allow-extra-chr --chr sex_linked_1
plink --bfile $PLINK --het --out outputs/sex_linked_2 --allow-extra-chr --chr sex_linked_2

# Window pi
vcftools --gzvcf $VCF --chr sex_linked_1 --keep males.txt \
         --window-pi 100000 --out sex_linked_1_males_pi

vcftools --gzvcf $VCF --chr sex_linked_1 --keep females.txt \
         --window-pi 100000 --out sex_linked_1_females_pi

vcftools --gzvcf $VCF --chr sex_linked_2 --keep males.txt \
         --window-pi 100000 --out sex_linked_2_males_pi

vcftools --gzvcf $VCF --chr sex_linked_2 --keep females.txt \
         --window-pi 100000 --out sex_linked_2_females_pi

# Compute missingnes
vcftools --gzvcf $VCF --chr outputs/sex_linked_1 --missing-indv --out sex_linked_1
vcftools --gzvcf $VCF --chr outputs/sex_linked_2 --missing-indv --out sex_linked_2