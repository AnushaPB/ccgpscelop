#!/bin/bash

source activate ccgpscelop

# INSTALLATION -------------------------------------------------------------------
# install smcpp if it doesn't exist already
if [ ! -d "smcpp" ]; then
    # sudo apt-get install -y python3-dev libgmp-dev libmpfr-dev libgsl0-dev
    # If the directory does not exist, execute the following commands
    git clone https://github.com/popgenmethods/smcpp.git
    cd smcpp
    pip install numpy scipy matplotlib pandas
    python setup.py install
    sudo apt install bedtools
    sudo apt install tabix
    cd ..
else
    echo "smcpp directory already exists."
fi

# CREATE UNCALLABLE SITES FILES -----------------------------------------------------------------------
# smcpp needs a bed file of sites that could not be called so that it properly identifies homozygous regions

# create paths
BASE_PATH="../../data"
PREFIX="58-Sceloporus"
#GENOME="$BASE_PATH/annotated_genome/jordan-uni4378-mb-hirise-65qvq__12-23-2023__final_assembly_relabelled.fasta.fai"
#GENOME_NAME="jordan-uni4378-mb-hirise-65qvq__12-23-2023__final_assembly_relabelled"
GENOME="$BASE_PATH/genome/rSceOcc1.20250428.p_ctg_chr_names.fna.gz"
GENOME_NAME="rSceOcc1.20250428.p_ctg_chr_names"
CALLABLE="$BASE_PATH/ccgp_data/58-Sceloporus_callable_sites.bed"

# Path to the uncompressed fasta file
GENOME_FA="${GENOME%.gz}"  # strip .gz extension

# If the .fai doesn't exist, generate it
if [ ! -f "${GENOME_FA}.fai" ]; then
    gunzip -c "$GENOME" > "$GENOME_FA"
    samtools faidx "$GENOME_FA"
fi

# create reverse masks from genome data if they don't exist
if [ -f "${PREFIX}_uncallable_sites.sorted.bed.gz" ]; then
    echo "Mask files already exists, skipping..."
else
    # Convert the .fai file to a BED file representing all genomic regions
    awk 'BEGIN {FS="\t"}; {print $1 "\t0\t" $2}' "${GENOME_FA}.fai" > "${GENOME_NAME}.bed"

    # Subtract the callable sites from the complete genomic regions
    bedtools subtract -a $GENOME_NAME.bed -b $CALLABLE > ${PREFIX}_uncallable_sites.bed

    # Index the bed files
    # Sort by chromosome (-k1,1) and start position (-k2,2n)
    # (required for proper indexing and downstream compatibility)
    sort -k1,1 -k2,2n ${PREFIX}_uncallable_sites.bed > ${PREFIX}_uncallable_sites.sorted.bed

    # Compress the sorted BED file using bgzip for tabix
    bgzip -c ${PREFIX}_uncallable_sites.sorted.bed > ${PREFIX}_uncallable_sites.sorted.bed.gz

    # Index the bed file
    tabix -p bed ${PREFIX}_uncallable_sites.sorted.bed.gz
fi

# FUNCTIONS FOR SMCPP -----------------------------------------------------------
# Function to process individual populations
process_population() {
    # Arguments:
    # pop: population name
    local pop=$1  

    # inds: Sample identifiers for the population
    local inds=$2    

    # Output directory path for storing the results of the population processing
    local output_dir="${BASE_PATH}/analysis/smcpp/${pop}"  

    # Creates the output directory if it does not exist
    mkdir -p "${output_dir}"

    # Pick a random individual from the population to distinguish the population
    # NOTE: the same distinguished individual must be used for all contigs in a population
    local distinguished_ind=$(echo $inds | tr ',' '\n' | shuf -n 1)

    # Write the distinguished individual to a file
    echo $distinguished_ind > "${output_dir}/${pop}_distinguished_ind.txt"

    # Loops through each contig in the contigs array
    # note: the ! operator is used for indirect reference, which means it gets the value of the variable whose name is stored in contigs
    for contig in "${CONTIGS[@]}"; do
        process_contig "$contig" "$output_dir" "$pop" "$inds" "$distinguished_ind"
    done

    # Select the *.smc.gz files listed in the output folder
    local smc_files=$(ls "${output_dir}"/*_"${pop}".smc.gz)

    # Estimates population parameters using smc++ estimate
    # and outputs the results to a specified file
    # --timepoints: This command specifies the starting and ending time points of the model. It accepts two numbers t1 tK specifying the starting and ending time points of the model (in generations). If not specified, SMC++ will use an heuristic to calculate the model time points points automatically.
    # timepoints are in generations, so we multiply by 2 assuming a generation time of 2 years
    # Sceloporus diverged ~1 Ma (1e6 years) ago = 1e6/2 = 5e5 generations
    # We will start at 100 generations ago (200 years ago)
    smc++ estimate -o "${output_dir}" "${MUTATION_RATE}" ${smc_files} --timepoints 1e2 5e5 > "${output_dir}/smcpp.out"

    # Generates a plot for the population analysis
    smc++ plot "${output_dir}/${pop}.pdf" "${output_dir}/model.final.json" -c
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
        smc++ vcf2smc -d ${distinguished_ind} ${distinguished_ind} \
                        "${VCF_FILE}" "${output_file}" "${contig}" \
                        "${pop}:${inds}" \
                        --mask "${MASK_FILE}"     
    fi
}

# RUNNING SMCPP -------------------------------------------------------------
# Define constants and file paths
BASE_PATH="../.."

# MUTATION RATE OPTION 1:
# mutation rate per site per year of 5.60 × 10 −10 , which was previously estimated for glass lizards (Anguidae Ophisaurus; Perry et al., 2018) .
# https://www.researchgate.net/publication/327964019_Molecular_Adaptations_for_Sensing_and_Securing_Prey_and_Insight_into_Amniote_Genome_Diversity_from_the_Garter_Snake_Genome
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
#VCF_FILE="${BASE_PATH}/data/ccgp_data/58-Sceloporus_complete_coords_annotated.vcf.gz"
VCF_FILE="${BASE_PATH}/data/ccgp_data/58-Sceloporus_complete_coords_annotated_occidentalis_only.vcf.gz"
MASK_FILE="${PREFIX}_uncallable_sites.sorted.bed.gz"

# Index the vcf file if it doesn't exist
if [ ! -f "${VCF_FILE}.csi" ]; then
    tabix -p vcf "${VCF_FILE}"
fi

# Make outputs file
mkdir -p outputs

# Get list of just chromosomes (not including sex-linked chromosomes)
CONTIG_FILE="outputs/chromosome_names.txt"
mapfile -t CONTIGS < "${CONTIG_FILE}"
# Print the contigs
echo ${CONTIGS[@]}

# Process each populationmkdir -p outputs
for j in {4..10}
do
  mkdir -p outputs/outputs$j

  for i in {1..8}
  do
    (
      inds=$(paste -sd, - < "${BASE_PATH}/analysis/admixture/outputs/k8_pop${i}.txt")
      process_population "pop${i}" "$inds" 2> pop${i}.stderr

      mv pop${i}/ outputs/outputs$j/pop${i}
    ) &
  done
  wait

  mv *stderr outputs/outputs$j
done

mv *bed* outputs

# Check the status of the background jobs
jobs
[1] 531019
[2] 531020
[3] 531021
[4] 531024
[5] 531026
[6] 531028
[7] 531032
[8] 531036