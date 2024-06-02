# Data processing of Sceloporus for ROH analyses
# Do the following in the analysis/anne directory

VCF="../../data/58-Sceloporus_maf05_minDP5_maxDP50_rmsamp60_mm80_rmsamp40.vcf.gz" # n=164 inds, not LD-pruned
FAI="../../data/genome/ncbi_dataset/data/GCA_023333645.1/GCA_023333645.1_rSceOcc1.0.p_genomic.fna.fai"

# ROH ------------------------------------------------------------------------------

# Define an array with your population identifiers
POPS_PATH="outputs"
POPS=("K5_pop1" "K5_pop2" "K5_pop3" "K5_pop4" "K5_pop5")

# Loop through each population
for POP in "${POPS[@]}"; do
    # Split by population
    bcftools view -S $POPS_PATH/${POP}.txt $VCF -Oz -o ${POP}.vcf.gz

    # Calculate ROH
    bcftools roh --threads 10 -G30 --AF-dflt 0.4 --estimate-AF GT,- -O z -o ${POP}.roh.gz ${POP}.vcf.gz
    zgrep RG ${POP}.roh.gz > ${POP}.rg.roh

    # Calculate fROH
    python -c "from roh import calc_roh; calc_roh('${POP}.rg.roh', '$FAI', '${POP}.froh')"

    # Move data into outputs directory
    mv *${POP}* outputs/
done