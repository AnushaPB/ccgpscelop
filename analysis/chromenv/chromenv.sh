#!/bin/bash
# Run in analysis/genetic_diversity
source activate ccgpscelop
BASE_PATH=../../data/ccgp_data
BFILE=$BASE_PATH/58-Sceloporus_complete_coords_annotated_chr

# For each chromosome 
for chr in {1..6}; do
  # Determine chromosome length 
  CHR_LENGTH=$(awk -v chr=$chr '$1==chr { if($4 > max) max=$4 } END { print max }' ${BFILE}.bim)

  # Make sure CHR_LENGTH and BFILE are set appropriately before running the loop
  for (( start=0; start<CHR_LENGTH; start+=500000 )); do
    end=$(( start + 500000 ))
    
    echo "Processing: chr $chr, from $start to $end"
    
    # Subset the data for the current bin and calculate heterozygosity
    plink --bfile $BFILE \
          --chr $chr \
          --from-bp $start \
          --to-bp $end \
          --het \
          --out outputs/output_chr${chr}_${start}_${end} \
          --allow-extra-chr
  done
done


# For each chromosome 
for chr in {1..6}; do
  # Determine chromosome length
  CHR_LENGTH=$(awk -v chr=$chr '$1==chr { if($4 > max) max=$4 } END { print max }' ${BFILE}.bim)
  
  for (( start=0; start<CHR_LENGTH; start+=500000 )); do
    end=$(( start + 500000 ))
    
    echo "Processing PCA: chr $chr, from $start to $end"
    
    # Subset the data for the current bin and run PCA (extracting first 2 PCs)
    plink --bfile $BFILE \
          --chr $chr \
          --from-bp $start \
          --to-bp $end \
          --pca 2 \
          --out outputs/pca_chr${chr}_${start}_${end} \
          --allow-extra-chr
  done
done