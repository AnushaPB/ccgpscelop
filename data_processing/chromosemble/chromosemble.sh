
# Install Satsuma
git clone https://github.com/bioinfologics/satsuma2.git
cd satsuma2
mkdir build
cd build
cmake ..
make
cd bin
# Add the bin directory 
export SATSUMA2_PATH=$(pwd)

# Move back into the chromosemble directory
cd ../../

# Define the paths to the genomes
OCC_GENOME=../../data/annotated_genome/jordan-uni4378-mb-hirise-65qvq__12-23-2023__final_assembly.fasta
UND_GENOME=../../data/genome_undulatus/ncbi_dataset/data/GCF_019175285.1/GCF_019175285.1_SceUnd_v1.1_genomic.fna

# Map scaffolds or contigs onto chromosome coordinates via synteny
$SATSUMA2_PATH/Chromosemble -t $UND_GENOME -q $OCC_GENOME -o outputs

# Align with nucmer
## Install mummer
wget https://github.com/mummer4/mummer/releases/download/v4.0.0rc1/mummer-4.0.0rc1.tar.gz
tar -xzf mummer-4.0.0rc1.tar.gz
rm mummer-4.0.0rc1.tar.gz
cd mummer-4.0.0rc1
./configure --prefix=$(pwd)
make
make install
cd ..
mummer-4.0.0rc1/nucmer  -p sceloporus  $UND_GENOME  $OCC_GENOME