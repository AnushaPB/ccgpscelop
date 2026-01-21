
# LANDFIRE Fire Regime Groups: https://www.landfire.gov/fire-regime/frg
wget -c -O .LF2016_FRG_200_CONUS.zip 'https://www.landfire.gov/data-downloads/US_200/LF2016_FRG_200_CONUS.zip'
unzip LF2016_FRG_200_CONUS.zip
rm ../../data/env/LF2016_FRG_200_CONUS.zip

# LANDFIRE Vegetation Departure: https://www.landfire.gov/vegetation/vdep
wget -c -O LF2016_VDep_200_CONUS.zip 'https://www.landfire.gov/data-downloads/US_200/LF2016_VDep_200_CONUS.zip'
unzip LF2016_VDep_200_CONUS.zip
rm LF2016_VDep_200_CONUS.zip

# LANDFIRE BPS (historical vegetation group)
wget -c -O LF2016_BPS_200_CONUS.zip 'https://www.landfire.gov/data-downloads/US_200/LF2016_BPS_200_CONUS.zip'
unzip LF2016_BPS_200_CONUS.zip
rm LF2016_BPS_200_CONUS.zip

# LANDFIRE Fire Return Interval
wget -c -O LF2016_FRI_200_CONUS.zip 'https://www.landfire.gov/data-downloads/US_200/LF2016_FRI_200_CONUS.zip'
unzip LF2016_FRI_200_CONUS.zip
rm LF2016_FRI_200_CONUS.zip

# LANDFIRE Percent Fire Severity
wget -c -O LF2016_PFS_200_CONUS.zip 'https://www.landfire.gov/data-downloads/US_200/LF2016_PFS_200_CONUS.zip'
unzip LF2016_PFS_200_CONUS.zip
rm LF2016_PFS_200_CONUS.zip

# LANDFIRE Historical Disturbance
wget -c -O LF2016_HDist_200_CONUS.zip 'https://www.landfire.gov/data-downloads/US_Disturbance/LF2016_HDist_200_CONUS.zip'
unzip LF2016_HDist_200_CONUS.zip
rm LF2016_HDist_200_CONUS.zip

# LANDFIRE current vegetation type
wget -c -O LF2016_EVT_200_CONUS.zip 'https://www.landfire.gov/data-downloads/US_200/LF2016_EVT_200_CONUS.zip'
unzip LF2016_EVT_200_CONUS.zip
rm LF2016_EVT_200_CONUS.zip