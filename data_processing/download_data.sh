#!/bin/bash

# List files in the specified directories
rsync --list-only hgdownload.soe.ucsc.edu::ccgp/58-Sceloporus_chr/QC/
rsync --list-only hgdownload.soe.ucsc.edu::ccgp/CCGP-module/58-Sceloporus_chr/
rsync --list-only hgdownload.soe.ucsc.edu::ccgp/CCGP-module/58-Sceloporus_chr/CCGP/

# Create directory for data
mkdir -p data/ccgp_data

# Download specific files to the data directory

# Callable sites:
rsync -avz hgdownload.soe.ucsc.edu::ccgp/58-Sceloporus_chr/58-Sceloporus_callable_sites.bed ../data/ccgp_data

# Coordinates:
rsync -avz hgdownload.soe.ucsc.edu::ccgp/58-Sceloporus_chr/QC/58-Sceloporus.coords.txt ../data/ccgp_data

# VCF files:
rsync -avz hgdownload.soe.ucsc.edu::ccgp/algatr/58-Sceloporus/58-Sceloporus_complete_coords_annotated.vcf.gz ../data/ccgp_data
rsync -avz hgdownload.soe.ucsc.edu::ccgp/algatr/58-Sceloporus/58-Sceloporus_complete_coords_annotated.vcf.gz.csi ../data/ccgp_data
rsync -avz hgdownload.soe.ucsc.edu::ccgp/algatr/58-Sceloporus/58-Sceloporus_complete_coords_pruned_0.6.vcf.gz ../data/ccgp_data
rsync -avz hgdownload.soe.ucsc.edu::ccgp/algatr/58-Sceloporus/58-Sceloporus_complete_coords_pruned_0.6.vcf.gz.csi ../data/ccgp_data

# Depth and missingness:
rsync -avz hgdownload.soe.ucsc.edu::ccgp/58-Sceloporus_chr/QC/58-Sceloporus.idepth ../data/ccgp_data
rsync -avz hgdownload.soe.ucsc.edu::ccgp/58-Sceloporus_chr/QC/58-Sceloporus.imiss ../data/ccgp_data

# RDA outputs:
mkdir -p ../analysis/gea/outputs
rsync -avz hgdownload.soe.ucsc.edu::ccgp/algatr/58-Sceloporus/RDA/ ../analysis/gea/outputs