import csv
from Bio import SeqIO

# Paths to files
gff_files = [
    "../../data/annotated_genome/annotation/cognate_trna.gff",
    "../../data/annotated_genome/annotation/non-coding.gff",
    "../../data/annotated_genome/annotation/protein-coding.w_func.gff"
]
fasta_files = [
    "../../data/annotated_genome/annotation/protein-coding.w_func.cds.fasta",
    "../../data/annotated_genome/annotation/protein-coding.w_func.prot.fasta",
    "../../data/annotated_genome/jordan-uni4378-mb-hirise-65qvq__12-23-2023__final_assembly.fasta"
]
mapping_file = "outputs/chromosome_labels.csv"

# Read the mapping file into a dictionary
def read_mapping(mapping_file):
    name_mapping = {}
    with open(mapping_file, "r") as map_file:
        reader = csv.reader(map_file)
        next(reader)  # Skip the header row
        for row in reader:
            old_name, new_name = row
            name_mapping[old_name] = new_name
    return name_mapping

# Update labels in GFF files
def update_gff_labels(input_file, output_file, name_mapping):
    with open(input_file, "r") as in_handle, open(output_file, "w") as out_handle:
        for line in in_handle:
            if line.startswith("#"):  # Skip comment lines
                out_handle.write(line)
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 9:
                # Update the sequence name (first column)
                seq_name = parts[0]
                if seq_name in name_mapping:
                    parts[0] = name_mapping[seq_name]
                out_handle.write('\t'.join(parts) + '\n')
            else:
                out_handle.write(line + '\n')

# Update labels in FASTA files
def update_fasta_labels(input_file, output_file, name_mapping):
    with open(input_file, "r") as input_handle, open(output_file, "w") as output_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            if record.id in name_mapping:
                record.id = name_mapping[record.id]
                #record.description = ""  # Optionally clear description if desired
            SeqIO.write(record, output_handle, "fasta")

def process_files():
    name_mapping = read_mapping(mapping_file)
    
    # Process GFF files
    for gff_file in gff_files:
        output_gff_file = gff_file.replace(".gff", "_relabelled.gff")
        update_gff_labels(gff_file, output_gff_file, name_mapping)
    
    # Process FASTA files
    for fasta_file in fasta_files:
        output_fasta_file = fasta_file.replace(".fasta", "_relabelled.fasta")
        update_fasta_labels(fasta_file, output_fasta_file, name_mapping)

if __name__ == "__main__":
    process_files()
