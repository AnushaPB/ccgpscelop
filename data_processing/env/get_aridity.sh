
# Download aridity index -----------------------------------------------------------------------

# Zomer, R.J.; Xu, J.; Trabuco, A. 2022. Version 3 of the Global Aridity Index and Potential Evapotranspiration Database. Scientific Data 9, 409. https://www.nature.com/articles/s41597-022-01493-1
# Data from figshare:
# https://doi.org/10.6084/m9.figshare.7504448 

curl -L \
  -A "Mozilla/5.0" \
  -o TEMP_figshare_download \
  "https://api.figshare.com/v2/file/download/56300327"

unzip TEMP_figshare_download

# Remove MACOSX hidden files
rm -rf __MACOSX

# Move the unzipped file to the data/env directory
mv Global-AI_ET0__annual_v3_1 ../data/env/Global-AI_ET0__annual_v3_1