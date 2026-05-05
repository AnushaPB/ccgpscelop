#!/bin/bash
# Example script to adapt the LASSI Plus workflow for one population across 11 chromosomes.
# Define the following environment variables (or modify the script directly):
#   VCF:    path to the input VCF file.
#   PREFIX: a prefix for the output files.
#   POPFILE: a text file with the sample IDs for the population of interest.
#
# For example, you can export these variables in your shell before running the script:
#   export VCF="mydata.vcf"
#   export PREFIX="mydata"
#   export POPFILE="pop.txt"

# Check that required variables are set
if [[ -z "$VCF" || -z "$PREFIX" || -z "$POPFILE" ]]; then
    echo "ERROR: Please set VCF, PREFIX, and POPFILE environment variables."
    exit 1
fi

########################################
# Step 1. Filter VCF for the chosen population
########################################
echo "Filtering VCF for population in $POPFILE..."
vcftools --vcf "$VCF" --keep "$POPFILE" --recode --out ${PREFIX}.pop
FILTERED_VCF="${PREFIX}.pop.recode.vcf"

########################################
# Step 2. Split the VCF into one VCF per chromosome (11 chromosomes)
########################################
echo "Separating VCF by chromosome..."
# Extract the header lines (all lines starting with "##")
grep "##" "$FILTERED_VCF" > ${PREFIX}_header.txt
# Extract the non-header lines (includes the column header)
grep -v "##" "$FILTERED_VCF" > ${PREFIX}_noheader.vcf

# Loop over chromosomes 1 to 11. (Adjust "scaffold_" prefix if needed.)
for chr in {1..11}; do
    echo "Processing chromosome scaffold_$chr..."
    awk -v chr="scaffold_$chr" '$1 == chr' ${PREFIX}_noheader.vcf > ${PREFIX}.chrom${chr}.vcf
done

########################################
# Step 3. Add the header back to each per-chromosome VCF
########################################
echo "Adding header to per-chromosome VCFs..."
# Extract the column header (first non-header line) from the filtered VCF.
grep -v "##" "$FILTERED_VCF" | awk 'NR==1{print}' > ${PREFIX}_column_header.txt
# Combine the meta-information header and column header.
cat ${PREFIX}_header.txt ${PREFIX}_column_header.txt > ${PREFIX}_main_header.vcf

# Prepend the header to each chromosome file.
for chr in {1..11}; do
    cat ${PREFIX}_main_header.vcf ${PREFIX}.chrom${chr}.vcf > ${PREFIX}.chrom${chr}.vcf.head
    mv ${PREFIX}.chrom${chr}.vcf.head ${PREFIX}.chrom${chr}.vcf
done

########################################
# Step 4. Run lassip on each chromosome
########################################
echo "Running lassip on each chromosome..."
for chr in {1..11}; do
    echo "Running lassip on chromosome scaffold_$chr..."
    lassip --vcf ${PREFIX}.chrom${chr}.vcf \
           --hapstats --winsize 201 --k 20 --calc-spec --winstep 100 \
           --out ${PREFIX}.chrom${chr} --pop "$POPFILE"
done

########################################
# Step 5. Compute average spectra to create a null background
########################################
echo "Computing average spectra for null background..."
lassip --spectra ${PREFIX}.chrom*.lassip.hap.spectra.gz --avg-spec --out ${PREFIX}.chromAll

########################################
# Step 6. Run saltiLASSI using the null spectra for each chromosome
########################################
echo "Running saltiLASSI on each chromosome..."
for chr in {1..11}; do
    echo "Running saltiLASSI on chromosome scaffold_$chr..."
    lassip --spectra ${PREFIX}.chrom${chr}.${PREFIX}.lassip.hap.spectra.gz \
           --salti --out ${PREFIX}.chrom${chr} \
           --null-spec ${PREFIX}.chromAll.lassip.null.spectra.gz
done

########################################
# Step 7. Remove column header from lassip output files and concatenate results
########################################
echo "Preparing final concatenated output..."
for chr in {1..11}; do
    awk 'NR>1' ${PREFIX}.chrom${chr}.lassip.hap.out > ${PREFIX}.chrom${chr}.lassip.hap.out.nohead
done

# Concatenate all output files (assumes file names sort in the correct order)
ls -v ${PREFIX}.chrom*.lassip.hap.out.nohead | xargs cat > ${PREFIX}_ALL.saltilassi.out

echo "All done. Final output saved to ${PREFIX}_ALL.saltilassi.out"
