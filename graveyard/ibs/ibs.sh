# Run in analysis/genetic_diversity
source activate ccgpscelop
BASE_PATH=../../data/ccgp_data_may2025backup
PLINK=$BASE_PATH/58-Sceloporus_complete_coords_annotated_chr

# Get long IBS tracts between individuals
plink --bfile $PLINK \
      --segment \
      --segment-snp 100 \
      --segment-kb 1000 \
      --out long_IBS_segments
