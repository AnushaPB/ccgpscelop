#!/bin/bash
# Run in data_processing/

# List files in the specified directories
rsync --list-only ccgp-download.gi.ucsc.edu::ccgp/58-Sceloporus/
rsync --list-only ccgp-download.gi.ucsc.edu::ccgp/CCGP-module/58-Sceloporus/

# Create directory for data
mkdir -p ../data/ccgp_data

# Download specific files to the data directory

# Callable sites:
rsync -avz ccgp-download.gi.ucsc.edu::ccgp/58-Sceloporus/58-Sceloporus_callable_sites.bed ../data/ccgp_data

# Coordinates:
rsync -avz ccgp-download.gi.ucsc.edu::ccgp/58-Sceloporus/QC/58-Sceloporus.coords.txt ../data/ccgp_data

# VCF files:
rsync -avz ccgp-download.gi.ucsc.edu::ccgp/CCGP-module/58-Sceloporus/58-Sceloporus_complete_coords_annotated.vcf.gz ../data/ccgp_data
rsync -avz ccgp-download.gi.ucsc.edu::ccgp/CCGP-module/58-Sceloporus/58-Sceloporus_complete_coords_annotated.vcf.gz.csi ../data/ccgp_data

rsync -avz ccgp-download.gi.ucsc.edu::ccgp/CCGP-module/58-Sceloporus/58-Sceloporus_annotated_pruned_0.6.vcf.gz ../data/ccgp_data
rsync -avz ccgp-download.gi.ucsc.edu::ccgp/CCGP-module/58-Sceloporus/58-Sceloporus_annotated_pruned_0.6.vcf.gz.csi ../data/ccgp_data

# Depth and missingness:
rsync -avz ccgp-download.gi.ucsc.edu::ccgp/58-Sceloporus/QC/58-Sceloporus.idepth ../data/ccgp_data
rsync -avz ccgp-download.gi.ucsc.edu::ccgp/58-Sceloporus/QC/58-Sceloporus.imiss ../data/ccgp_data

# RDA outputs:
mkdir -p ../analysis/gea/outputs
rsync --list-only ccgp-download.gi.ucsc.edu::ccgp/algatr/58-Sceloporus/RDA_no_pca/

rsync -avz ccgp-download.gi.ucsc.edu::ccgp/algatr/58-Sceloporus/RDA_no_pca/58-Sceloporus_scaledloadings.csv ../analysis/gea/outputs
rsync -avz ccgp-download.gi.ucsc.edu::ccgp/algatr/58-Sceloporus/RDA_no_pca/58-Sceloporus_unscaledloadings.csv ../analysis/gea/outputs
rsync -avz ccgp-download.gi.ucsc.edu::ccgp/algatr/58-Sceloporus/RDA_no_pca/58-Sceloporus_biplot.csv ../analysis/gea/outputs
rsync -avz ccgp-download.gi.ucsc.edu::ccgp/algatr/58-Sceloporus/RDA_no_pca/58-Sceloporus_RDA_outliers_full_Zscores.csv ../analysis/gea/outputs
rsync -avz ccgp-download.gi.ucsc.edu::ccgp/algatr/58-Sceloporus/RDA_no_pca/58-Sceloporus_RDA_cortest_full.csv  ../analysis/gea/outputs
