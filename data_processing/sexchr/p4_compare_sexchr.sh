# Confirmation that sex_linked_1 = SCAF_7 and sex_linked_2 = SCAF_8

# From Merly: https://drive.google.com/drive/folders/1bEsZPvrzxAkwjGJHOdhUchtIr7yTT5hs
ORIGINALGENOME=../../data/genome/rSceOcc1.20250428.p_ctg.fasta.gz
ANNOTATEDGENOME=../../data/genome/rSceOcc1.20250428.p_ctg_chr_names.fna.gz

seqkit fx2tab -n -l "$ORIGINALGENOME" \
  | awk 'BEGIN{OFS="\t"; print "scaffold","length_bp"} {print $1,$2}' \
  > outputs/rSceOcc1.20250428.p_ctg.scaffold_lengths.tsv
# SCAF_7	71284530
# SCAF_8	62430170

seqkit fx2tab -n -l "$ANNOTATEDGENOME" \
  | awk 'BEGIN{OFS="\t"; print "scaffold","length_bp"} {print $1,$2}' \
  > outputs/rSceOcc1.20250428.p_ctg_chr_names.scaffold_lengths.tsv
# sex_linked_1	71284530
# sex_linked_2	62430170

