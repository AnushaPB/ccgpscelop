BASE_PATH="../../data"
PREFIX="58-Sceloporus"
GENOME="$BASE_PATH/genome/ncbi_dataset/data/GCA_023333645.1/GCA_023333645.1_rSceOcc1.0.p_genomic.fna"
GENOME_NAME="GCA_023333645.1_rSceOcc1.0.p_genomic"
CALLABLE="$BASE_PATH/raw_data/58-Sceloporus_callable_sites.bed"

#  conda install -c bioconda genmap

#index reference genome
genmap index -F $GENOME -I index

#generate mappability files, give it an output directory 
# make more lenient k150
genmap map -K 150 -E 2 -I index -t -w -bg

#convert output into usable file for SNP filtering, this part is subjective re: how unique you want sites to be for mapping. 
awk '{for (i=$2+1; i<=$3 ; i++) print $1,i,$4}' DHP_genmap_k150_chrom_output.bedgraph > DHP.k150.chrom.genmap.tbl
awk '$3 > 0.5 {print $1"\t"$2}' DHP.k150.chrom.genmap.tbl > DHP.k150.chrom.genmap.pos