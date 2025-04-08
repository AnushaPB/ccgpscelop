# !/bin/bash
cd data/annotated_genome

input_file="jordan-uni4378-mb-hirise-65qvq__12-23-2023__final_assembly_relabelled.fasta.gz"
base_name="${input_file%.gz}"

# Decompress the file
gunzip -c "$input_file" > "$base_name"

# Create the index
samtools faidx "$base_name"