# envidatS3paths.txt was generated here: https://envicloud.wsl.ch/#/?prefix=chelsa%2Fchelsa_V2%2FGLOBAL%2F
# climatologies/1981-2010/bio/bio*
wget --no-host-directories --force-directories --input-file=envidatS3paths.txt
mkdir -p ../../data/env/chelsa 
mv chelsav2/GLOBAL/climatologies/1981-2010/bio/* ../../data/env/chelsa
rm -rf chelsav2