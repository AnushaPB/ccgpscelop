conda activate ccgpscelop

VCF=../data/ccgp_data/58-Sceloporus_complete_coords_annotated.vcf.gz
PRUNED=../data/ccgp_data/58-Sceloporus_annotated_pruned_0.6.vcf.gz
SAMPLEIDS=../data/final_sampleids.txt

# Retain final SampleIDs only
bcftools view -S ^$SAMPLEIDS $VCF -Oz -o ../data/ccgp_data/58-Sceloporus_complete_coords_annotated_final_samples.vcf.gz
bcftools view -S ^$SAMPLEIDS $PRUNED -Oz -o ../data/ccgp_data/58-Sceloporus_annotated_pruned_0.6_final_samples.vcf.gz

# Convert VCF to PLINK
plink --vcf ../data/ccgp_data/58-Sceloporus_complete_coords_annotated_final_samples.vcf.gz --make-bed --out ../data/ccgp_data/58-Sceloporus_complete_coords_annotated --allow-extra-chr --const-fid > vcf_to_plink.out 2> vcf_to_plink.err
plink --vcf ../data/ccgp_data/58-Sceloporus_annotated_pruned_0.6_final_samples.vcf.gz --make-bed --out ../data/ccgp_data/58-Sceloporus_annotated_pruned_0.6 --allow-extra-chr --const-fid > vcf_to_plink.out 2> vcf_to_plink.err

# Filter out chromosomes
# 70215227 out of 71154824 sites retained
# List scaffold names from PLINK file
awk '{print $1}' ../data/ccgp_data/58-Sceloporus_complete_coords_annotated.bim | sort | uniq > scaffold_names.txt
plink --bfile ../data/ccgp_data/58-Sceloporus_complete_coords_annotated --allow-extra-chr --chr 1 2 3 4 5 6 7 8 9 10 Scaffold_13__1_contigs__length_49873245 --make-bed --out ../data/ccgp_data/58-Sceloporus_complete_coords_annotated_chr

# full count: 50,181,443
# pruned count: 4,382,248
awk '{print $1}' ../data/ccgp_data/58-Sceloporus_annotated_pruned_0.6.bim | sort | uniq > scaffold_names_pruned.txt
plink --bfile ../data/ccgp_data/58-Sceloporus_annotated_pruned_0.6 --allow-extra-chr --chr 1 2 3 4 5 6 7 8 9 10 Scaffold_13__1_contigs__length_49873245 --make-bed --out ../data/ccgp_data/58-Sceloporus_annotated_pruned_0.6_chr