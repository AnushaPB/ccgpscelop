source activate cchpscelop

BASE_PATH=../../58-Sceloporus
VCF=$BASE_PATH/58-Sceloporus_maf05_minDP5_maxDP50_rmsamp60_mm80_rmsamp20.vcf.gz
GENOME="../../data/genome/ncbi_dataset/data/GCA_023333645.1/GCA_023333645.1_rSceOcc1.0.p_genomic.fna.fai"

# HETEROZYGOSITY ----------------------------------------------------------------------
# note: outputs homozygosity information
# FID: Family ID
# IID: Individual ID
# O(HOM): Observed number of homozygous genotypes
# O(HET): Observed number of heterozygous genotypes
# N(NM): Count of non-missing genotypes
# F: Inbreeding coefficient estimate
# calculate heterozygosity stats
plink --vcf $VCF --make-bed --out 58-Sceloporus --allow-extra-chr
plink --bfile 58-Sceloporus --het --out 58-Sceloporus --allow-extra-chr

# move files into outputs and delete unnecessary files
rm *.log
rm *.nosex
mv *.het outputs
mv *.bed outputs
mv *.bim outputs
mv *.fam outputs

# ROH -----------------------------
# Define an array with your population identifiers
POPS_PATH="../admixture/outputs"
#POPS=("pop1" "pop2" "pop3" "pop4" "pop5")
POPS="pop1"
# Loop through each population
for POP in "${POPS[@]}"; do
    # Split by population
    bcftools view -S $POPS_PATH/${POP}.txt $VCF -Oz -o ${POP}.vcf.gz

    # Calculate ROH
    bcftools roh --threads 10 -G30 --AF-dflt 0.4 --estimate-AF GT,- -O z -o ${POP}.roh.gz ${POP}.vcf.gz
    zgrep RG ${POP}.roh.gz > ${POP}.rg.roh

    # Calculate fROH
    python -c "from roh import calc_roh; calc_roh('${POP}.rg.roh', '\$GENOME', '${POP}.froh')"

    # Get maximum ROH length
    Rscript -e "print(max(read.table('${POP}.rg.roh')\$V6))"

    # Move data into output directory
    mkdir -p outputs
    mv *${POP}* outputs/
done
