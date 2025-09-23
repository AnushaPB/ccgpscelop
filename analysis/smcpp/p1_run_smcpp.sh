#!/bin/bash

# wrapper for calling smc++ inside Docker
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")"/../.. && pwd)"
VCF_FILE="${REPO_ROOT}/data/ccgp_data/58-Sceloporus_complete_coords_annotated_occidentalis_only.vcf.gz"
MASK_FILE="${REPO_ROOT}/analysis/smcpp/58-Sceloporus_uncallable_sites.sorted.bed.gz"
SMCPP_CMD="docker run --rm \
             -v ${REPO_ROOT}:${REPO_ROOT} \
             -w ${REPO_ROOT} \
             terhorst/smcpp:latest"

# CREATE UNCALLABLE SITES FILES -----------------------------------------------------------------------
# smcpp needs a bed file of sites that could not be called so that it properly identifies homozygous regions

# create paths
GENOME="${REPO_ROOT}/data/genome/rSceOcc1.20250428.p_ctg_chr_names.fna.gz"
CALLABLE="${REPO_ROOT}/data/ccgp_data/58-Sceloporus_callable_sites.bed"
GENOME_NAME="rSceOcc1.20250428.p_ctg_chr_names"

# Path to the uncompressed fasta file
GENOME_FA="${GENOME%.gz}"  # strip .gz extension

# If the .fai doesn't exist, generate it
if [ ! -f "${GENOME_FA}.fai" ]; then
    gunzip -c "$GENOME" > "$GENOME_FA"
    samtools faidx "$GENOME_FA"
fi

# create reverse masks from genome data if they don't exist
if [ -f "58-Sceloporus_uncallable_sites.sorted.bed.gz" ]; then
    echo "Mask files already exists, skipping..."
else
    # Convert the .fai file to a BED file representing all genomic regions
    awk 'BEGIN {FS="\t"}; {print $1 "\t0\t" $2}' "${GENOME_FA}.fai" > "${GENOME_NAME}.bed"

    # Subtract the callable sites from the complete genomic regions
    bedtools subtract -a $GENOME_NAME.bed -b $CALLABLE > 58-Sceloporus_uncallable_sites.bed

    # Index the bed files
    # Sort by chromosome (-k1,1) and start position (-k2,2n)
    # (required for proper indexing and downstream compatibility)
    sort -k1,1 -k2,2n 58-Sceloporus_uncallable_sites.bed > 58-Sceloporus_uncallable_sites.sorted.bed

    # Compress the sorted BED file using bgzip for tabix
    bgzip -c 58-Sceloporus_uncallable_sites.sorted.bed > 58-Sceloporus_uncallable_sites.sorted.bed.gz

    # Index the bed file
    tabix -p bed 58-Sceloporus_uncallable_sites.sorted.bed.gz
fi

# FUNCTIONS FOR SMCPP -----------------------------------------------------------

# Function to process individual populations
process_population() {
    #   pop:              population name (e.g. “pop1”)
    #   inds:             comma-separated list of all sample IDs in the population
    #   distinguished_ind: one sample ID from the provided list for this iteration
    #   it:               iteration index (1–10), used to namespace outputs

    local pop=$1
    local inds=$2
    local distinguished_ind=$3
    local it=$4

    # Output directory path for storing the results of the population processing
    # Lands in analysis/smcpp/it${it}/${pop}
    local output_dir="${REPO_ROOT}/analysis/smcpp/it${it}/${pop}"
    mkdir -p "${output_dir}"

    # Write the distinguished individual to a file
    echo "${distinguished_ind}" > "${output_dir}/${pop}_distinguished_ind.txt"

    # Loops through each contig in the CONTIGS array
    for contig in "${CONTIGS[@]}"; do
        process_contig "$contig" "$output_dir" "$pop" "$inds" "$distinguished_ind"
    done

    # Select the *.smc.gz files listed in the output folder
    smc_files=( "${output_dir}"/*_"${pop}".smc.gz )

   # Estimates population parameters using smc++ estimate
    # and outputs the results to a specified file
    # --timepoints: This command specifies the starting and ending time points of the model. It accepts two numbers t1 tK specifying the starting and ending time points of the model (in generations). If not specified, SMC++ will use an heuristic to calculate the model time points points automatically.
    # timepoints are in generations, so we multiply by 2 assuming a generation time of 2 years
    # Sceloporus diverged ~1 Ma (1e6 years) ago = 1e6/2 = 5e5 generations
    # We will start at 100 generations ago (200 years ago)
    $SMCPP_CMD estimate \
        -o "${output_dir}" \
        "${MUTATION_RATE}" \
        "${smc_files[@]}" \
        --timepoints 1e2 5e5 \
        > "${output_dir}/smcpp.out"

    # Generates a plot for the population analysis
    $SMCPP_CMD plot \
        "${output_dir}/${pop}.pdf" \
        "${output_dir}/model.final.json" \
        -c
}

# Function to process single contig
process_contig() {
    local contig=$1
    local output_dir=$2
    local pop=$3
    local inds=$4
    local distinguished_ind=$5

    # Defines the path for the output file for the current contig
    local output_file="${output_dir}/${contig}_${pop}.smc.gz"

    # If the output file doesn't exist, it generates a new one
    if [[ ! -f "${output_file}" ]]; then
        # Converts VCF to SMC format for the given contig and population
        $SMCPP_CMD vcf2smc -d ${distinguished_ind} ${distinguished_ind} \
                        "${VCF_FILE}" "${output_file}" "${contig}" \
                        "${pop}:${inds}" \
                        --mask "${MASK_FILE}"
    fi
}

# RUNNING SMCPP -------------------------------------------------------------

# MUTATION RATE OPTION 1:
# mutation rate per site per year of 5.60 × 10 −10 , which was previously estimated for glass lizards (Anguidae Ophisaurus; Perry et al., 2018) .
# https://www.researchgate.net/publication/327964019_Molecular_Adaptations_for_Sensing_and_Securing_Prey_and_Insight_into_Amniote_Genome_Diversity_from_THE_GARTER_SNAKE_Genome
# from paper: 5.60 × 10^−10 per site per year , generation time: 2-3 years
# converted to generation time: 5.60 × 10^−10 * 2 years = 11.2 x 10^-10 = 1.12e-9

# MUTATION RATE OPTION 2 (average across reptiles):
#  doi: 10.1038/s41586-023-05752-y
# "On average, mutation rates per generation are higher in reptiles (average of all species 1.17 × 10−8, 95% CI of the mean = 5.34 × 10−9 to 1.80 × 10−8)""
# AVERAGE MUTATION RATE: 1.17 × 10−8 per site per generation

# From Bouzid et al:
# Generation time of 2 years for S. occidentalis (Jameson & Allison, 1976).
# Bouzid et al. used a a human mutation rate of 1.0 × 10−8 substitutions per site per generation
MUTATION_RATE=1e-8

# Files
VCF_FILE="${REPO_ROOT}/data/ccgp_data/58-Sceloporus_complete_coords_annotated_occidentalis_only.vcf.gz"
MASK_FILE="${REPO_ROOT}/analysis/smcpp/58-Sceloporus_uncallable_sites.sorted.bed.gz"
CALLABLE="${REPO_ROOT}/data/ccgp_data/58-Sceloporus_callable_sites.bed"

# Index the vcf file if it doesn't exist
if [ ! -f "${VCF_FILE}.csi" ]; then
    tabix -p vcf "${VCF_FILE}"
fi

# Get list of just chromosomes (not including sex-linked chromosomes)
CONTIG_FILE="outputs/chromosome_names.txt"
mapfile -t CONTIGS < "${CONTIG_FILE}"
# Print the contigs
echo ${CONTIGS[@]}

# Make output directory
mkdir -p outputs

# Define variables
NUM_ITS=10
ADMIX_DIR="${REPO_ROOT}/analysis/admixture/outputs"

# iterate 1…10 (“it”)
for it in $(seq 1 $NUM_ITS); do
  for i in {1..8}; do
    (
      pop="pop${i}"

      # inds: full comma-separated list for this population
      inds=$(paste -sd, - < "${ADMIX_DIR}/k8_${pop}.txt")

      # pick the it-th distinguished individual 
      dist_file="${ADMIX_DIR}/k8_${pop}_distinguished_individuals.txt"
      distinguished_ind=$(sed -n "${it}p" "${dist_file}")
      echo "Processing population: $pop, distinguished individual: $distinguished_ind, iteration: $it"

      process_population "$pop" "$inds" "$distinguished_ind" "$it" 2> "${pop}.stderr"
    ) &
  done

  wait
  mv "it${it}" outputs/
  mv *.stderr outputs/it${it}/
done

[1] 469632
[2] 469633
[3] 469634
[4] 469636
[5] 469639
[6] 469642
[7] 469645
[8] 469648