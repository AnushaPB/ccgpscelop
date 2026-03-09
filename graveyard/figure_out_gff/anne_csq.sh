conda activate ccgpscelop

# Install agat to handle GFF files if not already installed
conda install -c bioconda agat
GENOME=../../data/genome/rSceOcc1.20250428.p_ctg_chr_names.fna
ANNOTATION=../../data/genome/annotation/complete.genomic.gff
VCF=../../data/ccgp_data/58-Sceloporus_complete_coords_annotated_occidentalis_only.vcf.gz

# First, try with bcftools gff2gff script
git clone https://github.com/samtools/bcftools.git
gzip -c $ANNOTATION > complete.genomic.gff.gz
zcat complete.genomic.gff.gz | ./bcftools/misc/gff2gff | bgzip -c > complete.genomic.fixed.gff.gz
bcftools csq -f "$GENOME" -g complete.genomic.fixed.gff.gz --local-csq "$VCF" -Ob -o csq.bcf

# Now, try with agat
# Standardize the NCBI gff by converting to gtf
agat_convert_sp_gff2gtf.pl --gff $ANNOTATION -o ncbi.tmp.gtf
# Convert back to gff3
agat_convert_sp_gtf2gff3.pl --gtf ncbi.tmp.gtf -o ncbi.normalized.gff

# Mapping file was downloaded from Ensembl BioMart
# Params to click: old_id > new_id
# PodMur_mart_export.txt is mapped using Podarcis muralis genome; Chicken_mart_export.txt is using chicken
# Convert from NCBI to ENSEMBL
agat_sp_manage_attributes.pl --gff ncbi.normalized.gff --map PodMur_mart_export.txt --attribute ID -o ensembl.tmp.gff
agat_sp_manage_attributes.pl --gff ensembl.tmp.gff --map PodMur_mart_export.txt --attribute Parent -o ensembl.ids_fixed.gff

# Without true Ensembl IDs
agat_convert_sp_gff2gtf.pl -gff ncbi.gff -o tmp.gtf
agat_convert_sp_gtf2gff3.pl -gtf tmp.gtf -o ncbi.cleaned.gff

# Extract all fields (no filtering) to a file
# bcftools query -f '%CHROM\t%POS\t%INFO/BCSQ\n' outputs/csq.bcf > outputs/all_variant_csq.txt
# wc -l outputs/all_variant_csq.txt

# Remove the temporary file
# rm outputs/fixed_annotation.gff