

awk '$1=="chr6" {sum += $3 - $2} END {print sum}' 58-Sceloporus_callable_sites.bed 
# 113360742

awk '$1=="chr6"' 58-Sceloporus_callable_sites.bed | head -n 1
# chr6    519234  519270

tabix 58-Sceloporus_complete_coords_annotated.vcf.gz chr6 | head -n 1 | awk '{print $1, $2}'
# chr6 24519693


# CHECK IN data/notannotated_data/raw_data
# JALMGF010000011.1 aligns with the first 25 Mb of chr6


awk '$1=="JALMGF010000011.1" {sum += $3 - $2} END {print sum}' 58-Sceloporus_callable_sites.bed 
# 22952338

tabix 58-Sceloporus_filteredQC.vcf.gz JALMGF010000011.1 | head -n 1 | awk '{print $1, $2}'
# 
# (No JALMGF010000011.1)

# Confirmed correct code with:
tabix 58-Sceloporus_filteredQC.vcf.gz JALMGF010000001.1 | head -n 1 | awk '{print $1, $2}'
# 10000001.1 44461JALMGF010000001.1 7578

bcftools index 58-Sceloporus_filteredQC.vcf.gz  
bcftools index 58-Sceloporus.pruned.vcf.gz

# Get list of scaffolds in genome that aren't in the vcf
GENOME=data/notannotated_data/genome/ncbi_dataset/data/GCA_023333645.1/GCA_023333645.1_rSceOcc1.0.p_genomic.fna
VCF=data/notannotated_data/raw_data_ccgp/QC/58-Sceloporus_filteredQC.vcf.gz

grep "^>" "$GENOME" | sed 's/^>//' | cut -d' ' -f1 > genome_scaffolds.txt
bcftools view -h "$VCF" | grep "^##contig" | sed 's/##contig=<ID=//' | cut -d',' -f1 > vcf_scaffolds.txt

comm -23 <(sort genome_scaffolds.txt) <(sort vcf_scaffolds.txt) > nonoverlap.txt

rsync -avz hgdownload.soe.ucsc.edu::ccgp/58-Sceloporus/58-Sceloporus_raw.vcf.gz .
rsync -avz hgdownload.soe.ucsc.edu::ccgp/58-Sceloporus/58-Sceloporus_raw.vcf.gz.tbi .


# Extract chromosome lengths from VCF header
bcftools view -h $VCF | grep '##contig=' | sed -E 's/.*ID=([^,]+),length=([0-9]+).*/\1\t\2/' > scaffold_lengths.txt

bcftools view $VCF -Ou | bcftools query -f '%CHROM\t%POS\t%POS\t%DP\n'  | awk -v OFS='\t' '{{print $1,$2-1,$2,$4}}' > mean_depth.txt
