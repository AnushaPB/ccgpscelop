# CALCULATE DATA SIZES

rsync -avzn --stats --list-only hgdownload.soe.ucsc.edu::ccgp/CCGP-module
# convert
bytes=6493995835684
terabytes=$(echo "scale=2; $bytes / 1024^4" | bc)
echo $terabytes

# Create list of species folders
rsync -avzn --include="*/" --exclude="*" --list-only hgdownload.soe.ucsc.edu::ccgp/CCGP-module/ | \
# Extract the species folders
grep '/CCGP' | \
awk '{print $5}' | \
# Remove the /CCGP from the end of the folder name
sed 's|/CCGP||' | \
# Exclude the weird folders
grep -v -e 'cali_small' -e 'no_filter' -e 'old' -e 'original' -e '.bak' -e 'extended' > CCGP-module_species

for species in `cat CCGP-module_species`
do 
  rsync -avzn --stats --list-only hgdownload.soe.ucsc.edu::ccgp/CCGP-module/$species/${species}_clean_snps.vcf.gz
done

for species in $(cat CCGP-module_species)
do 
  size=$(rsync -avzn --stats --list-only hgdownload.soe.ucsc.edu::ccgp/CCGP-module/$species/${species}_clean_snps.vcf.gz | awk '/Total file size/ {print $4}')
  size=$(echo $size | tr -d ',')  # remove commas
  echo "$species,$size" >> species_sizes.csv
done

# read in the sizes and sum to get size in terms of terabytes
python3 -c "
import pandas as pd
sizes = pd.read_csv('species_sizes.csv', names=['species', 'size'])
total_size_tb = round(sizes['size'].sum() / (1024 ** 4), 2)
nprojects = sizes[sizes['size'] > 0].shape[0]
print(f'Total size of all files: {total_size_tb} TB')
print(f'Number of projects: {nprojects}')
"

# DOWNLOAD DATA
# Create list of species folders
rsync -avzn --include="*/" --exclude="*" --list-only hgdownload.soe.ucsc.edu::ccgp/CCGP-module/ | \

# Extract the species folders
grep '/CCGP' | \
awk '{print $5}' | \
# Remove the /CCGP from the end of the folder name
sed 's|/CCGP||' | \
# Exclude the weird folders
grep -v -e 'cali_small' -e 'no_filter' -e 'old' -e 'original' -e '.bak' -e 'extended' -e 'plink' > CCGP-module_species

mkdir -p clean_snps
# Download the clean snps files
for species in `cat CCGP-module_species`
do 
  rsync -avz hgdownload.soe.ucsc.edu::ccgp/CCGP-module/$species/${species}_clean_snps.vcf.gz clean_snps
done

mkdir -p coords
# Download the coordinates
for species in `cat CCGP-module_species`
do 
  rsync -avz hgdownload.soe.ucsc.edu::ccgp/$species/QC/*coords.txt coords
done

