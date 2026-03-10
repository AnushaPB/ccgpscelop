# Run in analysis/genetic_diversity
source activate ccgpscelop
BASE_PATH=../../data/ccgp_data
# !!! NOTE: running into issues with this file so currently het/window_pi was created on Rancor using Anne's 58-Sceloporus_annotated plink files
PLINK=$BASE_PATH/58-Sceloporus_complete_coords_annotated
PLINK_CHR=$BASE_PATH/58-Sceloporus_complete_coords_annotated_chr
PLINK_PRUNED=$BASE_PATH/58-Sceloporus_annotated_pruned_0.6_chr

# 1) Quick subset of chr1 with light QC
plink --bfile $PLINK_CHR --chr 1 --make-bed --out chr1.quick

# 2) Randomly keep ~X% of SNPs to cut pair count way down
plink --bfile chr1.quick --thin 0.50 --make-bed --out chr1.quick.thin

# 3) Fast LD: 500 kb window, keep only pairs with r² ≥ 0.2 (much fewer rows)
#    Add --threads N if you want
plink --bfile chr1.quick.thin \
  --r2 gz --ld-window-kb 500 --ld-window 999999 --ld-window-r2 0.2 \
  --out chr1.fastld


PLINK_PATH=../../analysis/genetic_diversity/chr1.quick
# ./bin/LDBlockShow \
#   -InPlink $PLINK_PATH \
#   -OutPut chr1_LD \
#   -Region 1:1000000-2000000 \
#   -SeleVar 2 \
#   -OutPdf

plink --bfile $PLINK_PATH --indep-pairwise 1000 100 0.1 --out chr1_pruned
plink --bfile $PLINK_PATH --extract chr1_pruned.prune.in --make-bed --out chr1_thin
./bin/LDBlockShow -InPlink chr1_thin -OutPut chr1_LD -SeleVar 2 -OutPdf
END=$(awk '$1=="1"{if($4>m)m=$4}END{print m+0}' chr1_thin.bim)


mkdir -p out
CHR=1   # or chr1
END=$(awk -v C="$CHR" '$1==C{if($4>m)m=$4}END{print m+0}' chr1_thin.bim)
#END=1000000
./bin/LDBlockShow \
  -InPlink chr1_thin \
  -OutPut out/chr1_LD_full \
  -Region ${CHR}:1-${END} \
  -SeleVar 2 \
  -BlockType 5 \
  -OutPng

# # (optional) also make a PDF from the SVG without touching ImageMagick policy:
# sudo apt-get install -y librsvg2-bin
# rsvg-convert -f pdf -o out/chr1_LD_full.pdf out/chr1_LD_full.svg


# NEW:
plink --bfile $PLINK_PATH --indep-pairwise 10000 100 0.8 --out chr1_pruned_loose
plink --bfile $PLINK_PATH --extract chr1_pruned_loose.prune.in --make-bed --out chr1_loose

END=$(awk -v C="$CHR" '$1==C{if($4>m)m=$4}END{print m+0}' chr1_thin.bim)
END=100000
./bin/LDBlockShow -InPlink chr1_loose -OutPut out/chr1_loose_full -Region 1:1-$END -SeleVar 2 -BlockType 5 -OutPng
