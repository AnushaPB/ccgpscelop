conda activate ccgpscelop

# Run to use GEE to get NDVI/gHM data
# You need a Google Earth Engine account to do this and you need to set up your credentials with the gee package
python get_ndvi.py
python get_ghm.py
# Download NDVI/gHM data from google drive and put in data/env
https://drive.google.com/drive/u/0/folders/1bs8guawRXfEQvdriDPzMMZ8lFaBk0pC7

# Download prism data
bash get_prism.sh

# Download and process chelsa data
bash get_chelsa.sh
Rscript process_chelsa.R