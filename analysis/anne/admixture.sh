# Data processing of Sceloporus for ADMIXTURE
# Do the following in the analysis/anne directory

PREFIX="58-Sceloporus_maf05_minDP5_maxDP50_rmsamp60_mm80_rmsamp40_r60"
VCF="../../data/58-Sceloporus_maf05_minDP5_maxDP50_rmsamp60_mm80_rmsamp40_r60.vcf.gz" # n=164 inds, LD-pruned

# ADMIXTURE ------------------------------------------------------------------------------

# Create folder for outputs
mkdir -p outputs

# Make plink bed, bim, and fam files required for admixture
plink --vcf $VCF --out $PREFIX --allow-extra-chr --make-bed --const-fid

# Fix chromosome numbers
awk '{$1=0;print $0}' $PREFIX.bim > $PREFIX.bim.tmp # make temp file
mv $PREFIX.bim.tmp $PREFIX.bim # replace original file with temp file

# Run ADMIXTURE
for K in 1 2 3 4 5 6 7 8 9 10; do admixture --cv $PREFIX.bed $K | tee $PREFIX.${K}.out; done
mv *58-Sceloporus* outputs

# Quickly retrieve CV scores for output
grep -h CV outputs/*.out

# Two individuals will be removed subsequently: single individual in Mojave: IW1426 (Scelocci_IW1426); 
# individual that has seemingly incorrect coords: Scelocci_BUR76

# Run admixture.Rmd which will generate files with individuals assigned to each cluster and
# save them as e.g., K5_pop1.txt, K5_pop2.txt, etc.

