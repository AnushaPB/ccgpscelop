#!/bin/bash

source activate ccgpscelop

# INSTALLATION -------------------------------------------------------------------
# install smcpp if it doesn't exist already
if [ ! -d "smcpp" ]; then
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
GENOME="$BASE_PATH/genome/ncbi_dataset/data/GCA_023333645.1/GCA_023333645.1_rSceOcc1.0.p_genomic.fna.fai"
GENOME_NAME="GCA_023333645.1_rSceOcc1.0.p_genomic"
CALLABLE="$BASE_PATH/raw_data/58-Sceloporus_callable_sites.bed"

# create reverse masks from genome data if they don't exist
if [ -f "${PREFIX}_uncallable_sites.sorted.bed.gz" ]; then
    echo "Mask files already exists, skipping..."
else
    # Convert the .fai file to a BED file representing all genomic regions
    awk 'BEGIN {FS="\t"}; {print $1 "\t0\t" $2}' $GENOME > $GENOME_NAME.bed

    # Subtract the callable sites from the complete genomic regions
    bedtools subtract -a $GENOME_NAME.bed -b $CALLABLE > ${PREFIX}_uncallable_sites.bed

    # Index the bed files
    sort -k1,1 -k2,2n ${PREFIX}_uncallable_sites.bed > ${PREFIX}_uncallable_sites.sorted.bed
    bgzip -c ${PREFIX}_uncallable_sites.sorted.bed > ${PREFIX}_uncallable_sites.sorted.bed.gz
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

    # Loops through each contig in the  contigs array
    # note: the ! operator is used for indirect reference, which means it gets the value of the variable whose name is stored in contigs
    for contig in "${CONTIGS[@]}"; do
        process_contig "$contig" "$output_dir" "$pop" "$inds"
    done

    # SUBSET 10 CONTIGS: Select the first 10 *.smc.gz files listed in the output folder
    local smc_files=$(ls "${output_dir}"/*_"${pop}".smc.gz | head -n 10)

    # Estimates population parameters using smc++ estimate
    # and outputs the results to a specified file
    smc++ estimate -o "${output_dir}" "${MUTATION_RATE}" "${output_dir}"/*_"${pop}".smc.gz --timepoints 1e2 1e5 > "${output_dir}/smcpp.out"
    #smc++ estimate -o "${output_dir}" "${MUTATION_RATE}" "${output_dir}"/*_"${pop}".smc.gz > "${output_dir}/smcpp.out"

    # Generates a plot for the population analysis
    smc++ plot "${output_dir}/${pop}.pdf" "${output_dir}/model.final.json" -c
}

# Function to process single contig
process_contig() {
    local contig=$1
    local output_dir=$2
    local pop=$3
    local inds=$4

    # Defines the path for the output file for the current contig
    local output_file="${output_dir}/${contig}_${pop}.smc.gz"

    # If the output file doesn't exist, it generates a new one
    if [[ ! -f "${output_file}" ]]; then
        # Converts VCF to SMC format for the given contig and population
        smc++ vcf2smc "${VCF_FILE}" "${output_file}" "${contig}" \
                          "${pop}:${inds}" \
                          --mask "${MASK_FILE}"
    fi
}

# RUNNING SMCPP -------------------------------------------------------------
# Define constants and file paths
BASE_PATH="../.."

# mutation rate per site per year of 5.60 × 10 −10 , which was previously estimated for glass lizards (Anguidae Ophisaurus; Perry et al., 2018) .
# https://www.researchgate.net/publication/327964019_Molecular_Adaptations_for_Sensing_and_Securing_Prey_and_Insight_into_Amniote_Genome_Diversity_from_the_Garter_Snake_Genome
# from paper: 5.60 × 10^−10 per site per year , generation time: 2-3 years
# converted to generation time: 5.60 × 10^−10 * 2 years = 11.2 x 10^-10 = 1.12e-9
MUTATION_RATE=1.12e-9

# Files
VCF_FILE="${BASE_PATH}/data/processed_data/58-Sceloporus_maf05_minDP5_maxDP50_rmsamp60_mm80_rmsamp40.vcf.gz"
MASK_FILE="${PREFIX}_uncallable_sites.sorted.bed.gz"

# Index the vcf file
tabix -p vcf "${VCF_FILE}"

# Extract contig information for contigs with length > 1Mb
CONTIG_FILE="outputs/contigs.txt"
if [ -f $CONTIF_FILE ]; then
    mapfile -t CONTIGS_ALL < "${CONTIG_FILE}"
else
    mapfile -t CONTIGS_ALL < <(zgrep '##contig=' "${VCF_FILE}" | awk -F'[<>,]' '{for(i=1;i<=NF;i++) {if($i ~ /ID=/) id=substr($i,4); if($i ~ /length=/) len=substr($i,8);} if(len+0 > 1000000) print id}')
    # write contigs_all out to eht CONTIG FILE
    printf "%s\n" "${CONTIGS_ALL[@]}" > "${CONTIG_FILE}"
fi

# Subset contigs
#CONTIGS=("${CONTIGS_ALL[0]}")
CONTIGS=("${CONTIGS_ALL[@]:0:10}")
# Print the contigs
echo ${CONTIGS[@]}

# Process each population
for i in {1..7}
do
  eval inds_pop$i=$(cat "${BASE_PATH}/analysis/admixture/outputs/k7_pop$i.txt" | paste -sd, -)
  pop="inds_pop${i}"
  process_population "pop${i}" "${!pop}" 2> pop$i.stderr &
done
wait

# Check the status of the background jobs
jobs
# Move files into outputs
mv *bed* outputs
mv *stderr outputs
for i in {1..7}
do
  mv pop$i/ outputs/pop$i
done