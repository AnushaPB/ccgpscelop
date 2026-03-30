# Download genomic data
bash p0_download_genome.sh

# Download environmental data
bash env/JOB_ENV.sh

# Index genome
bash p1_index_genome.sh

# Write text file with final sample IDs
Rscript p2_write_sample_ids.R

# Make final VCFs and convert VCF to PLINK format
bash p3_vcf_to_plink.sh

# Thin data
bash p4_thin_data.sh

# Count callable sites
bash p5_count_callable.sh

# Convert GFF format for csq
Rscript p6_gff_conversion.R

# Run csq to get consequences of SNPs
bash p7_csq.sh

# Convert csq files to bed files
Rscript p8_csq_to_bed.R

# NOTE: 
# chromsemble directory was used to annotate scaffolds with chromosome numbers
# sexchr directory was used to identify sex chromosome scaffolds
