source activate ccgpscelop
# conda install -c bioconda admixture

# run from analysis/admixture

# Create folder for outputs
mkdir -p outputs

# Specify prefix for file formatting
PREFIX="58-Sceloporus_pruned_0.6_thinned_10kb_chr"
# based on code from: https://github.com/ccgproject/ccgpWorkflow/blob/227e6506cf03274be7295e9de6411bb34162d97c/workflow/modules/qc/Snakefile#L126
# Copy plink files over from processed_data
cp ../../data/ccgp_data/$PREFIX.bed .
cp ../../data/ccgp_data/$PREFIX.bim .
cp ../../data/ccgp_data/$PREFIX.fam .

# Remove beckii individuals
awk '{print "0", $1}' ../../data/beckii_sampleids.txt > ../../data/beckii_sampleids_plink.txt
plink --bfile $PREFIX --remove ../../data/beckii_sampleids_plink.txt --make-bed --out $PREFIX.tmp
mv $PREFIX.tmp.bed $PREFIX.bed
mv $PREFIX.tmp.bim $PREFIX.bim
mv $PREFIX.tmp.fam $PREFIX.fam

# Fix chromosome numbers
awk '{$1=0;print $0}' $PREFIX.bim > $PREFIX.bim.tmp # make temp file
mv $PREFIX.bim.tmp $PREFIX.bim # replace original file with temp file

# Run ADMIXTURE
for K in {1..15}; do
  admixture --cv $PREFIX.bed $K | tee $PREFIX.${K}.out
done
mv *$PREFIX* outputs/