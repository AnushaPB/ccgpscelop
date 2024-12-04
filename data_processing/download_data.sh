#!/bin/bash

rsync --list-only hgdownload.soe.ucsc.edu::ccgp/58-Sceloporus_chr/QC/
rsync --list-only hgdownload.soe.ucsc.edu::ccgp/CCGP-module/58-Sceloporus_chr/
rsync --list-only hgdownload.soe.ucsc.edu::ccgp/CCGP-module/58-Sceloporus_chr/CCGP/

mkdir -p data/ccgp_data

rsync -avz hgdownload.soe.ucsc.edu::ccgp/58-Sceloporus_chr/58-Sceloporus_callable_sites.bed data/ccgp_data

rsync -avz hgdownload.soe.ucsc.edu::ccgp/58-Sceloporus_chr/QC/58-Sceloporus.coords.txt data/ccgp_data

rsync -avz hgdownload.soe.ucsc.edu::ccgp/CCGP-module/58-Sceloporus_chr/58-Sceloporus_clean_snps.vcf.gz data/ccgp_data
#rsync --list-only hgdownload.soe.ucsc.edu::ccgp/CCGP-module/58-Sceloporus_chr/58-Sceloporus_clean_snps.vcf.gz

rsync -avz hgdownload.soe.ucsc.edu::ccgp/CCGP-module/58-Sceloporus_chr/CCGP/58-Sceloporus_annotated_pruned_0.6.vcf.gz data/ccgp_data
#rsync --list-only hgdownload.soe.ucsc.edu::ccgp/CCGP-module/58-Sceloporus_chr/CCGP/58-Sceloporus_annotated_pruned_0.6.vcf.gz

rsync -avz hgdownload.soe.ucsc.edu::ccgp/CCGP-module/58-Sceloporus_chr/CCGP/58-Sceloporus_annotated_pruned_0.6.dist data/ccgp_data
rsync -avz hgdownload.soe.ucsc.edu::ccgp/CCGP-module/58-Sceloporus_chr/CCGP/58-Sceloporus_annotated_pruned_0.6.dist.id data/ccgp_data

rsync -avz hgdownload.soe.ucsc.edu::ccgp/58-Sceloporus_chr/QC/58-Sceloporus.idepth data/ccgp_data
rsync -avz hgdownload.soe.ucsc.edu::ccgp/58-Sceloporus_chr/QC/58-Sceloporus.imiss data/ccgp_data

rsync -avz hgdownload.soe.ucsc.edu::ccgp/CCGP-module/58-Sceloporus_chr/CCGP/58-Sceloporus_pruned_mil.vcf.gz data/ccgp_data

rsync -avz hgdownload.soe.ucsc.edu::ccgp/58-Sceloporus_chr/QC/58-Sceloporus.het data/ccgp_data