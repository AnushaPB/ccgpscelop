# Download CHELSA current bioclimatic variables -----------------------------------------------------------------------
# chelsa_climatologies.txt was generated here: https://envicloud.wsl.ch/#/?prefix=chelsa%2Fchelsa_V2%2FGLOBAL%2F
# climatologies/1981-2010/bio/bio*
wget --no-host-directories --force-directories --input-file=chelsa_climatologies.txt
mkdir -p ../../data/env/chelsa 
mv chelsav2/GLOBAL/climatologies/1981-2010/bio/* ../../data/env/chelsa
rm -rf chelsav2

# Download PMIP3 data -----------------------------------------------------------------------------------------------
# Downscaled global climatological data from the last glacial maximum (21.000 years ago). The CHELSA LGM data is based on a implementation of the CHELSA algorithm on PMIP3 data. Currently we provide data for seven PMIP3 GCMs. NCAR-CCSM4, MRI-CGCM3, MPI-ESM-P, MIROC-ESM, CESS-FGOALS-g2, IPSL-CM5A-LR, CNRM-CM5.
# The Community Climate System Model (CCSM) is a coupled general circulation model (GCM) developed by the University Corporation for Atmospheric Research (UCAR) with funding from the National Science Foundation (NSF), the Department of Energy (DoE), and the National Aeronautics and Space Administration (NASA).
wget --no-host-directories --force-directories --input-file=chelsa_pmip.txt
mkdir -p ../../data/env/chelsa_pmip
mv chelsav1/pmip3/bioclim/* ../../data/env/chelsa_pmip
rm -rf chelsav1

# Download TraCE21k data ----------------------------------------------------------------------------------------------
# CHELSA-TraCE21k data provides monthly climate data for temperature and precipitation at 30 arcsec spatial resolution in 100-yeartime steps for the last 21,000 years. Paleo orography at high spatial resolution and at each timestep is created by combining high resolution information on glacial cover from current and Last Glacial Maximum (LGM) glacier databases with the interpolation of a dynamic ice sheet model (ICE6G) and a coupling to mean annual temperatures from CCSM3-TraCE21k. Based on the reconstructed paleo orography, mean annual temperature and precipitation was downscaled using the CHELSA V1.2 algorithm.
wget --no-host-directories --force-directories --input-file=chelsa_trace21k.txt
mkdir -p ../../data/env/chelsa_trace21k
mv chelsav1/chelsa_trace/bio/* ../../data/env/chelsa_trace21k
rm -rf chelsav1

# Download CHELSA future bioclimatic variables -----------------------------------------------------------------------
wget --no-host-directories --force-directories --input-file=chelsa_future.txt
mkdir -p ../../data/env/chelsa_future
mv chelsav2/GLOBAL/climatologies/2071-2100/GFDL-ESM4/ssp585/bio/* ../../data/env/chelsa_future
#rm -rf chelsav2

# Download CHELSA vpd -----------------------------------------------------------------------
wget --no-host-directories --force-directories --input-file=chelsa_vpd.txt
mkdir -p ../../data/env/chelsa_vpd
mv chelsa02/chelsa/global/climatologies/vpd/1981-2010/* ../../data/env/chelsa_vpd
rm -rf chelsa02

wget --no-host-directories --force-directories --input-file=chelsa_cmi.txt
mkdir -p ../../data/env/chelsa_cmi
mv chelsa02/chelsa/global/climatologies/cmi/1981-2010/* ../../data/env/chelsa_cmi
rm -rf chelsa02