# Download historical fire regime data from LANDFIRE

# Download contemporary fire history metrics for the conterminous US from 1984-2020 from USGS
# Link from (will expire but can be regenerated) : https://www.sciencebase.gov/catalog/item/6244bbeed34e21f8276030bf
wget "https://prod-is-usgs-sb-prod-content.s3.amazonaws.com/6244bbeed34e21f8276030bf/conus_1984_2020_metrics.zip?AWSAccessKeyId=AKIAI7K4IX6D4QLARINA&Expires=1767223870&Signature=fY2NwZc1bIQYcSVya5I4HL3mkno%3D" \
-O conus_1984_2020_metrics.zip
unzip conus_1984_2020_metrics.zip -d conus_1984_2020_metrics
mv conus_1984_2020_metrics ../data/env/conus_fire_history_metrics_1984_2020
rm conus_1984_2020_metrics.zip

# LANDFIRE Fire Regime Groups: https://www.landfire.gov/fire-regime/frg
wget -c -O .LF2016_FRG_200_CONUS.zip 'https://www.landfire.gov/data-downloads/US_200/LF2016_FRG_200_CONUS.zip'
unzip LF2016_FRG_200_CONUS.zip
rm LF2016_FRG_200_CONUS.zip
mv LF2016_FRG_200_CONUS ../data/env/LF2016_FRG_200_CONUS

# LANDFIRE Vegetation Departure: https://www.landfire.gov/vegetation/vdep
wget -c -O LF2016_VDep_200_CONUS.zip 'https://www.landfire.gov/data-downloads/US_200/LF2016_VDep_200_CONUS.zip'
unzip LF2016_VDep_200_CONUS.zip
rm LF2016_VDep_200_CONUS.zip
mv LF2016_VDep_200_CONUS ../data/env/LF2016_VDep_200_CONUS

# LANDFIRE BPS (historical vegetation group)
wget -c -O LF2016_BPS_200_CONUS.zip 'https://www.landfire.gov/data-downloads/US_200/LF2016_BPS_200_CONUS.zip'
unzip LF2016_BPS_200_CONUS.zip
rm LF2016_BPS_200_CONUS.zip
mv LF2016_BPS_200_CONUS ../data/env/LF2016_BPS_200_CONUS

# LANDFIRE Fire Return Interval
wget -c -O LF2016_FRI_200_CONUS.zip 'https://www.landfire.gov/data-downloads/US_200/LF2016_FRI_200_CONUS.zip'
unzip LF2016_FRI_200_CONUS.zip
rm LF2016_FRI_200_CONUS.zip
mv LF2016_FRI_200_CONUS ../data/env/LF2016_FRI_200_CONUS

# LANDFIRE Percent Fire Severity
wget -c -O LF2016_PFS_200_CONUS.zip 'https://www.landfire.gov/data-downloads/US_200/LF2016_PFS_200_CONUS.zip'
unzip LF2016_PFS_200_CONUS.zip
rm LF2016_PFS_200_CONUS.zip
mv LF2016_PFS_200_CONUS ../data/env/LF2016_PFS_200_CONUS

# LANDFIRE Historical Disturbance
wget -c -O LF2016_HDist_200_CONUS.zip 'https://www.landfire.gov/data-downloads/US_Disturbance/LF2016_HDist_200_CONUS.zip'
unzip LF2016_HDist_200_CONUS.zip
rm LF2016_HDist_200_CONUS.zip
mv LF2016_HDist_200_CONUS ../data/env/LF2016_HDist_200_CONUS

# LANDFIRE current vegetation type
wget -c -O LF2016_EVT_200_CONUS.zip 'https://www.landfire.gov/data-downloads/US_200/LF2016_EVT_200_CONUS.zip'
unzip LF2016_EVT_200_CONUS.zip
rm LF2016_EVT_200_CONUS.zip
mv LF2016_EVT_200_CONUS ../data/env/LF2016_EVT_200_CONUS