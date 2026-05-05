# Run in analysis/genetic_diversity
source activate ccgpscelop
BASE_PATH=../../data/ccgp_data
# !!! NOTE: running into issues with this file so currently het/window_pi was created on Rancor using Anne's 58-Sceloporus_annotated plink files
PLINK=$BASE_PATH/58-Sceloporus_complete_coords_annotated
PLINK_PRUNED=$BASE_PATH/58-Sceloporus_annotated_pruned_0.6

sudo apt-get update
sudo apt-get install -y build-essential git zlib1g-dev

git clone --recurse-submodules https://github.com/williamslab/ibis.git
cd ibis
make


./ibis \
  -bfile $PLINK \   # PLINK prefix (no .bed/.bim/.fam extension)
  -noConvert \             # treat Mb as cM (no genetic map present)
  -min_l 2 \               # min tract length ≥2 Mb
  -mt 500 \                # ≥500 markers per tract
  -er 0.002 \              # allow ~1 mismatch per tract (0.002 × 500 ≈ 1)
  -threads 8 \             # use 8 CPU threads
  -o ibis_result           # output prefix: ibis_, etc.