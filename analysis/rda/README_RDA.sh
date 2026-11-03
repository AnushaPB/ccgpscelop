# The following provides information to run the RDA analysis, which was run using a snakemake pipeline

# Code takes in the following files:
#       data/ccgp_data/58-Sceloporus_annotated_pruned_0.6.vcf.gz which is the LD-pruned vcf
#       data/genome/rSceOcc1.20250428.p_ctg_chr_names.fna.fai which is reference genome file
#       

#   Snakefile:  To see the order of operations
#   config.yaml: To see parameters altered for Sceloporus
#   Scripts can be found in analysis/rda/algatr
#       gather_rda.R
#       general_functions.R pulls in env and gen data, processes, etc.
#       generate_contig_files.py detects and splits vcf by contig
#       RDA.R code to actually run RDA on each scaffold
#       setup_algatr_env.R sets up algatr env (installs all required packages for algatr)

# RDA outputs all saved in analysis/gea/outputs/RDA_results/

# =============== Running the snakemake pipeline ===============
# Activate environment
conda activate snakemake

# Run the snakemake pipeline
snakemake --dir ../../ --profile profiles/adaptive -n