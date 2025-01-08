#!/bin/bash
wget -r ftp://prism.oregonstate.edu/monthly/tmean/

mkdir -p TEMP
mv prism.oregonstate.edu/monthly/tmean/* TEMP/

# Unzip all .zip files within each folder, keeping the directory structure
for year in {1895..2024}
do
  find "TEMP/$year" -type f -name "*.zip" -exec unzip -o -d "TEMP/$year" {} \;
done

# Move all files to the data directory
mkdir -p ../../data/env/prism/PRISM_tmean_stable_4kmM3
mv TEMP/* ../../data/env/prism/PRISM_tmean_stable_4kmM3

# Remove the TEMP directory and the prism.oregonstate.edu directory
rm -r TEMP prism.oregonstate.edu
