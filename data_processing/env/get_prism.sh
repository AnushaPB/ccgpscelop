#!/bin/bash

# Source: https://prism.oregonstate.edu/historical/

# Download historical (1895-2025) monthly tmean temperature from PRISM
# Note: this downloads both annual and monthly mean data
# Coding of name:
# time_series = time series data
# us = United States
# an = an” refers to “all networks” and is the default PRISM data. The other option is “lt” which refers to “long term,”
# which is an 800m monthly version of PRISM data that incorporates observations from only long-term
# established weather stations. 
# 800m = 800m resolution
# tmean = mean temperature
wget -r ftp://prism.oregonstate.edu/time_series/us/an/800m/tmean/monthly

# Make temporary directory
mkdir -p TEMP

# Move data into temporary directory
mv prism.oregonstate.edu/time_series/us/an/800m/tmean/monthly/* TEMP/

# Unzip all .zip files within each folder, keeping the directory structure
for year in {1895..2024}
do
  find "TEMP/$year" -type f -name "*.zip" -exec unzip -o -d "TEMP/$year" {} \;
done

# Move all files to the data directory
mkdir -p ../../data/env/prism/prism_time_series_us_an_800m_tmean
mv TEMP/* ../../data/env/prism/prism_time_series_us_an_800m_tmean

# Remove the TEMP directory and the prism.oregonstate.edu directory
rm -r TEMP prism.oregonstate.edu
