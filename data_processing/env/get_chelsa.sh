# Download CHELSA current bioclimatic variables
# chelsa_climatologies.txt was generated here: https://envicloud.wsl.ch/#/?prefix=chelsa%2Fchelsa_V2%2FGLOBAL%2F
# climatologies/1981-2010/bio/bio*
wget --no-host-directories --force-directories --input-file=chelsa_climatologies.txt
mkdir -p ../../data/env/chelsa 
mv chelsav2/GLOBAL/climatologies/1981-2010/bio/* ../../data/env/chelsa
rm -rf chelsav2

# Download PMIP data
wget --no-host-directories --force-directories --input-file=chelsa_pmip.txt
mkdir -p ../../data/env/chelsa_pmip
mv chelsav1/pmip3/bioclim/* ../../data/env/chelsa_pmip
rm -rf chelsav1

# Download TraCE21k data
wget --no-host-directories --force-directories --input-file=chelsa_trace21k.txt
mkdir -p ../../data/env/chelsa_trace21k
mv chelsav1/chelsa_trace/bio/* ../../data/env/chelsa_trace21k
rm -rf chelsav1