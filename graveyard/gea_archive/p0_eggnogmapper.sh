#!/usr/bin/env bash

### 2) Create & activate a new env with eggNOG-mapper ###
conda activate ccgpscelop
conda install eggnog-mapper diamond hmmer sqlite samtools seqkit

### 3) Download the eggNOG data (~20 GB) ###
download_eggnog_data.py  --data_dir outputs/eggnog -y

### 4) Extract your seven genes from the master FASTA ###
# (make sure you run this from your project root)
FASTA=data/genome/annotation/complete.proteins.faa
OUT=subset_genes.faa

# a) list of your locus_tags
cat > genes_of_interest.txt <<EOF
egapxtmp_013352
egapxtmp_009456
egapxtmp_020386
egapxtmp_013358
egapxtmp_013367
egapxtmp_013373
egapxtmp_000232
EOF

# b) index & extract with samtools faidx
samtools faidx "${FASTA}"
samtools faidx "${FASTA}" genes_of_interest.txt > "${OUT}"

# (Alternatively, using seqkit if your headers aren’t simple names)
# seqkit grep -n -f genes_of_interest.txt "${FASTA}" > "${OUT}"

### 5) Run eggNOG-mapper on your subset ###
emapper.py \
  -i "${OUT}" \
  --output emapper_results/subset_emapper \
  --data_dir /path/to/eggnog_data \
  --cpu 8 \
  --go_evidence non-electronic \
  --override

echo "Done!  Results in emapper_results/subset_emapper.emapper.annotations"
