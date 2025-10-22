import sys

# Convert GFF to a format bcftools csq can parse (Ensembl-style IDs and structure)
def convert_gff(input_file, output_file):
    gene_map = {}        # maps original gene IDs -> new gene:ID
    transcript_map = {}  # maps original transcript IDs -> new transcript:ID
    
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        outfile.write("##gff-version 3\n")  # header required by bcftools csq
        
        for line in infile:
            if line.startswith("#") or not line.strip():
                continue  # skip comments and blank lines
            
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue  # ignore malformed lines
            
            chrom, source, feature_type, start, end, score, strand, phase, attributes = fields
            # turn "ID=xxx;Parent=yyy" into a dict
            attr_dict = {kv.split("=")[0]: kv.split("=")[1] for kv in attributes.split(";") if "=" in kv}
            
            if feature_type == "gene":
                gene_id = attr_dict.get("ID", "unknown_gene")
                new_gene_id = f"gene:{gene_id}"
                gene_map[gene_id] = new_gene_id
                # output with prefixed ID and biotype
                outfile.write("\t".join([chrom, source, "gene", start, end, score, strand, phase,
                                         f"ID={new_gene_id};biotype=protein_coding\n"]))
            
            elif feature_type == "mRNA":
                transcript_id = attr_dict.get("ID", "unknown_transcript")
                parent_gene = attr_dict.get("Parent", "unknown_gene")
                new_transcript_id = f"transcript:{transcript_id}"
                transcript_map[transcript_id] = new_transcript_id
                # look up parent gene with new prefix
                new_parent = gene_map.get(parent_gene, f"gene:{parent_gene}")
                outfile.write("\t".join([chrom, source, "transcript", start, end, score, strand, phase,
                                         f"ID={new_transcript_id};Parent={new_parent};biotype=protein_coding\n"]))
            
            elif feature_type in ["exon", "CDS", "three_prime_UTR", "five_prime_UTR"]:
                parent_transcript = attr_dict.get("Parent", "unknown_transcript")
                # point to transcript:ID instead of raw Parent
                new_parent = transcript_map.get(parent_transcript, f"transcript:{parent_transcript}")
                new_attributes = f"Parent={new_parent}"
                if feature_type == "CDS":
                    new_attributes += f";Phase={phase}"  # CDS requires phase
                outfile.write("\t".join([chrom, source, feature_type, start, end, score, strand, phase,
                                         new_attributes + "\n"]))

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python convert_gff.py input.gff output.gff")
    else:
        convert_gff(sys.argv[1], sys.argv[2])
        print(f"Conversion complete. Output saved to {sys.argv[2]}")
