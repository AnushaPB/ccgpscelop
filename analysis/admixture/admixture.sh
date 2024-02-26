source activate ccgpscelop
# conda install -c bioconda admixture

# run from analysis/admixture

# Create folder for outputs
mkdir -p outputs

# VCF file
VCF=../../58-Sceloporus/58-Sceloporus_maf05_minDP5_maxDP50_rmsamp60_mm80_rmsamp20_r60.vcf.gz

# based on code from: https://github.com/ccgproject/ccgpWorkflow/blob/227e6506cf03274be7295e9de6411bb34162d97c/workflow/modules/qc/Snakefile#L126
# make plink files
plink --vcf $VCF --out 58-Sceloporus --allow-extra-chr --make-bed --const-fid --noweb

# fix chromosome numbers
awk '{$1=0;print $0}' 58-Sceloporus.bim > 58-Sceloporus.bim.tmp # make temp file
mv 58-Sceloporus.bim.tmp 58-Sceloporus.bim # replace original file with temp file
# run admixture
for i in {2..10}
do
 admixture --cv 58-Sceloporus.bed $i > 58-Sceloporus.log${i}.out
done
# move data into subfolder
mv *58-Sceloporus* outputs
