# HETEROZYGOSITY ----------------------------------------------------------------------
# note: outputs homozygosity information
# FID: Family ID
# IID: Individual ID
# O(HOM): Observed number of homozygous genotypes
# O(HET): Observed number of heterozygous genotypes
# N(NM): Count of non-missing genotypes
# F: Inbreeding coefficient estimate
# calculate heterozygosity stats
# use const-fid flag to set FID (population ID) to 0; otherwise Error: Multiple instances of '_' in sample ID.
for VCF in ../clean_snps/*_clean_snps.vcf.gz
do
  # Extract the base name of the file to use in output files
  BASENAME=$(basename "$VCF" _clean_snps.vcf.gz)

  # calculate heterozygosity stats
  plink --vcf $VCF --make-bed --out $BASENAME --allow-extra-chr --const-fid
  plink --bfile $BASENAME --het --out $BASENAME --allow-extra-chr
done

# note: calculated the same as https://www.cell.com/current-biology/pdfExtended/S0960-9822(21)00548-0: 
# "we define heterozygosity as the fraction of heterozygous genotypes out of all
# called genotypes (including invariant sites) within an individual"

# Calculate the number of sites in the bed file
for BED in ../bed/*.bed
do
  BASENAME=$(basename "$BED" _callable_sites.bed)
  NSITES=$(wc -l $BED | awk '{print $1}')
  echo "$BASENAME,$NSITES" >> bed_nsites.csv
done