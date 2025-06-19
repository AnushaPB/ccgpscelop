conda activate ccgpscelop
# Install Satsuma
git clone https://github.com/bioinfologics/satsuma2.git
cd satsuma2
mkdir build
cd build
cmake ..
make
cd bin

# Move back into the chromosemble directory
cd ../../../

export SATSUMA2_PATH=satsuma2/build/bin

# Gunzip the genome files
gunzip ../../data/new_genome/rSceOcc1.20250428.a_ctg.fasta.gz
gunzip ../../data/new_genome/rSceOcc1.20250428.p_ctg.fasta.gz

# Define the paths to the genomes
OCC_GENOME_A=../../data/new_genome/rSceOcc1.20250428.a_ctg.fasta
OCC_GENOME_P=../../data/new_genome/rSceOcc1.20250428.p_ctg.fasta
OLD_OCC_GENOME=../../data/annotated_genome/jordan-uni4378-mb-hirise-65qvq__12-23-2023__final_assembly_relabelled.fasta
UND_GENOME=../../data/genome_undulatus/ncbi_dataset/data/GCF_019175285.1/GCF_019175285.1_SceUnd_v1.1_genomic.fna

# Map scaffolds or contigs onto chromosome coordinates via synteny
$SATSUMA2_PATH/Chromosemble -t $UND_GENOME -q $OCC_GENOME_A -o outputs_a
$SATSUMA2_PATH/Chromosemble -t $UND_GENOME -q $OCC_GENOME_P -o outputs_p
$SATSUMA2_PATH/Chromosemble -t $OCC_GENOME_P -q $OLD_OCC_GENOME -o new_vs_old

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

# Subset the genomes to the first 12 largest scaffolds
conda install -c bioconda seqkit
seqkit fx2tab -l -n $UND_GENOME | sort -k2,2nr | head -n 12 | cut -f1 > top12_und.txt
seqkit fx2tab -l -n $OCC_GENOME_P | sort -k2,2nr | head -n 12 | cut -f1 > top12_occ.txt
seqkit grep -n -f top12_und.txt $UND_GENOME > und_top12.fasta
seqkit grep -n -f top12_occ.txt $OCC_GENOME_P > occ_top12.fasta

mummer-4.0.0rc1/nucmer -p sceloporus -t 12 und_top12.fasta occ_top12.fasta
delta-filter -1 sceloporus.delta > sceloporus.filter
mummerplot --fat --layout --filter --size large --png -p sceloporus sceloporus.filter



# MUMMER
# Get top 20 scaffolds >1Mb from UND genome
seqkit fx2tab -l -n "$OLD_OCC_GENOME" | awk '$2 >= 1000000' | sort -k2,2nr | head -n 20 | cut -f1 > top20_occ_old.txt

# Get top 20 scaffolds >1Mb from OCC genome
seqkit fx2tab -l -n "$OCC_GENOME_P" | awk '$2 >= 1000000' | sort -k2,2nr | head -n 20 | cut -f1 > top20_occ_p.txt

# Subset each genome to top 20 scaffolds
seqkit grep -n -f top20_occ_old.txt "$OLD_OCC_GENOME" > top20_occ_old.fasta
seqkit grep -n -f top20_occ_p.txt "$OCC_GENOME_P" > top20_occ_p.fasta

# Run nucmer with 20 threads
mummer-4.0.0rc1/nucmer -p sceloporus -t 20 top20_occ_old.fasta top20_occ_p.fasta

# Filter for best 1-to-1 alignments
delta-filter -1 sceloporus.delta > sceloporus.filter

# Plot alignments
mummerplot --fat --layout --filter --size large --png -p sceloporus sceloporus.filter
