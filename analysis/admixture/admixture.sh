source activate ccgpscelop
# conda install -c bioconda admixture

# run from analysis/admixture

# Create folder for outputs
mkdir -p outputs

# Specify prefix for file formatting
PREFIX="../../58-Sceloporus/58-Sceloporus_maf05_minDP5_maxDP50_rmsamp60_mm80_rmsamp20_r60"

# based on code from: https://github.com/ccgproject/ccgpWorkflow/blob/227e6506cf03274be7295e9de6411bb34162d97c/workflow/modules/qc/Snakefile#L126
# Make plink files
plink --vcf $PREFIX.vcf.gz --out $PREFIX --allow-extra-chr --make-bed --const-fid --noweb

# Fix chromosome numbers
awk '{$1=0;print $0}' $PREFIX.bim > $PREFIX.bim.tmp # make temp file
mv $PREFIX.bim.tmp $PREFIX.bim # replace original file with temp file

# Do two separate runs of ADMIXTURE
mkdir -p outputs/admixture_k1-10
for K in 1 2 3 4 5 6 7 8 9 10; do admixture --cv $PREFIX.bed $K | tee $PREFIX.${K}.out; done
mv *58-Sceloporus* outputs/admixture_k1-10

mkdir -p outputs/admixture_k3-15
for K in 3 4 5 6 7 8 9 10 11 12 13 14 15; do admixture --cv $PREFIX.bed $K | tee $PREFIX.${K}.out; done
mv *58-Sceloporus* outputs/admixture_k3-15