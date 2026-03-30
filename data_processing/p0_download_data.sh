#!/bin/bash
# Run in data_processing/

# List files in the specified directories
rsync --list-only ccgp-download.gi.ucsc.edu::ccgp/58-Sceloporus/
rsync --list-only ccgp-download.gi.ucsc.edu::ccgp/CCGP-module/58-Sceloporus/

# Create directory for data
mkdir -p ../data/ccgp_data

# Download specific files to the data directory

# LIst files
rsync --list-only ccgp-download.gi.ucsc.edu::ccgp/CCGP-module/58-Sceloporus/

# Callable sites:
rsync -avz ccgp-download.gi.ucsc.edu::ccgp/58-Sceloporus/58-Sceloporus_callable_sites.bed ../data/ccgp_data

# Coordinates:
rsync -avz ccgp-download.gi.ucsc.edu::ccgp/58-Sceloporus/QC/58-Sceloporus.coords.txt ../data/ccgp_data

# VCF files:
rsync -avz ccgp-download.gi.ucsc.edu::ccgp/CCGP-module/58-Sceloporus/58-Sceloporus_complete_coords_annotated.vcf.gz ../data/ccgp_data
rsync -avz ccgp-download.gi.ucsc.edu::ccgp/CCGP-module/58-Sceloporus/58-Sceloporus_complete_coords_annotated.vcf.gz.csi ../data/ccgp_data

rsync -avz ccgp-download.gi.ucsc.edu::ccgp/CCGP-module/58-Sceloporus/58-Sceloporus_annotated_pruned_0.6.vcf.gz ../data/ccgp_data
rsync -avz ccgp-download.gi.ucsc.edu::ccgp/CCGP-module/58-Sceloporus/58-Sceloporus_annotated_pruned_0.6.vcf.gz.csi ../data/ccgp_data

rsync -avz ccgp-download.gi.ucsc.edu::ccgp/CCGP-module/58-Sceloporus/58-Sceloporus_clean_snps.vcf.gz ../data/ccgp_data
rsync -avz ccgp-download.gi.ucsc.edu::ccgp/CCGP-module/58-Sceloporus/58-Sceloporus_clean_snps.vcf.gz.tbi ../data/ccgp_data

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

# Other data sources:
# Genome:
# https://drive.google.com/drive/folders/1bEsZPvrzxAkwjGJHOdhUchtIr7yTT5hs

# Annotated genome:
# https://drive.google.com/drive/folders/1wHUzbCK9Ep-cCRcqVcGhYKIjfeIZ008q?usp=share_link

# NCBI genome:
mkdir -p genome
cd genome
curl -L -o genome.zip "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCA_023333645.1/download?include_annotation_type=GENOME_FASTA"
unzip genome.zip
rm genome.zip
# index genome
samtools faidx ncbi_dataset/data/GCA_023333645.1/GCA_023333645.1_rSceOcc1.0.p_genomic.fna 
cd ..

# NCBI undulatus genome:
mkdir -p genome_undulatus
cd genome_undulatus
curl -L -o genome_undulatus.zip "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_019175285.1/download?include_annotation_type=GENOME_FASTA"
unzip genome_undulatus.zip
rm genome_undulatus.zip
samtools faidx ncbi_dataset/data/GCF_019175285.1/GCF_019175285.1_SceUnd_v1.1_genomic.fna