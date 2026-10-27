
#!/bin/bash

# Run from scripts/data_processing/
# Copy files for Dryad archive

DATA=../data
DRYAD="$DATA/dryad"
CCGP="$DATA/ccgp_data"
GENOME="$DATA/genome"

mkdir -p "$DRYAD"

# VCF
cp "$CCGP/58-Sceloporus_complete_coords_annotated_final_samples.vcf.gz" "$DRYAD/"

# PLINK datasets (BED, BIM, FAM)
cp "$CCGP/58-Sceloporus_complete_coords_annotated_chr."{bed,bim,fam} "$DRYAD/"
cp "$CCGP/58-Sceloporus_annotated_pruned_0.6_chr."{bed,bim,fam} "$DRYAD/"
cp "$CCGP/58-Sceloporus_pruned_0.6_thinned_10kb_chr."{bed,bim,fam} "$DRYAD/"

# Sampling coordinates and callable sites
cp "$CCGP/58-Sceloporus.coords.txt" "$DRYAD/"
cp "$CCGP/58-Sceloporus_callable_sites.bed" "$DRYAD/"

# Reference genome and annotation
cp "$GENOME/rSceOcc1.20250428.p_ctg_chr_names.fna" "$DRYAD/"
cp "$GENOME/rSceOcc1.20250428.p_ctg_chr_names.fna.fai" "$DRYAD/"
cp "$GENOME/annotation/complete.genomic.gff" "$DRYAD/"

# Copy README file
cp "README_dryad.md" "$DRYAD/README.md"

# Zip the Dryad directory
cd "$DATA"
zip -r dryad.zip dryad
