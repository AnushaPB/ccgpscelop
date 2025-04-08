#!/bin/bash
# Run in analysis/genetic_diversity
source activate ccgpscelop
BASE_PATH=../../data/ccgp_data
BFILE=$BASE_PATH/58-Sceloporus_complete_coords_annotated_chr

# NOTE THIS FILE IS GIGANTIC
plink --bfile $BFILE --chr 5 --r2 --ld-window-kb 100 --ld-window 100 --ld-window-r2 0 --out ld_chr5 --allow-extra-chr
