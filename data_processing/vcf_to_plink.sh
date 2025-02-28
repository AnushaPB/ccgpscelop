conda activate ccgpscelop

VCF=../data/ccgp_data/58-Sceloporus_complete_coords_annotated.vcf.gz
PRUNED=../data/ccgp_data/58-Sceloporus_complete_coords_pruned_0.6.vcf.gz

# Convert VCF to PLINK
plink --vcf $VCF --make-bed --out ../data/ccgp_data/58-Sceloporus_complete_coords_annotated --allow-extra-chr --const-fid > vcf_to_plink.out 2> vcf_to_plink.err

plink --vcf $PRUNED --make-bed --out ../data/ccgp_data/58-Sceloporus_complete_coords_pruned_0.6 --allow-extra-chr --const-fid > vcf_to_plink.out 2> vcf_to_plink.err

#plink --vcf 58-Sceloporus_complete_coords_annotated.vcf.gz --make-bed --out 58-Sceloporus_complete_coords_annotated --allow-extra-chr --const-fid

#plink --vcf 58-Sceloporus_complete_coords_pruned_0.6.vcf.gz --make-bed --out 58-Sceloporus_complete_coords_pruned_0.6 --allow-extra-chr --const-fid

# Filter out chromosomes
# 70215227 out of 71154824 
# List scaffold names from PLINK file
awk '{print $1}' ../data/ccgp_data/58-Sceloporus_complete_coords_annotated.bim | sort | uniq > scaffold_names.txt

plink --bfile ../data/ccgp_data/58-Sceloporus_complete_coords_annotated --allow-extra-chr --chr 1 2 3 4 5 6 7 8 9 10 11 Scaffold_13__1_contigs__length_49873245 --make-bed --out ../data/ccgp_data/58-Sceloporus_complete_coords_annotated_chr

# Before count: 50,181,443
# After count: 4,382,248
plink --bfile ../data/ccgp_data/58-Sceloporus_complete_coords_pruned_0.6 --allow-extra-chr --chr 1 2 3 4 5 6 7 8 9 10 11 Scaffold_13__1_contigs__length_49873245 --make-bed --out ../data/ccgp_data/58-Sceloporus_complete_coords_pruned_0.6_chr