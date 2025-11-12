# Run in analysis/genetic_diversity
source activate ccgpscelop
BASE_PATH=../../data/ccgp_data
PLINK_CHR=$BASE_PATH/58-Sceloporus_complete_coords_annotated_chr
OUT_DIR="outputs"

# Loop over all bin files (TSV, one SNP ID per line)
for f in "${OUT_DIR}"/variants_MAF_*.tsv; do
  [ -e "$f" ] || { echo "No variants_MAF_*.tsv files found in ${EXTRACT_DIR}"; break; }

  # Make a safe tag for filenames (replace anything non-alnum with underscores)
  base="$(basename "$f" .tsv)"
  tag="$(printf '%s' "$base" | sed 's/[^A-Za-z0-9._-]/_/g')"   # e.g., variants_MAF_0.00_0.01

  echo ">>> Running PLINK --het for bin list: $f"
  plink \
    --bfile "${PLINK_CHR}" \
    --extract "$f" \
    --het \
    --allow-extra-chr \
    --out "${OUT_DIR}/het_${tag}"
done


NONSYN=outputs/allnonsyn

for f in "${OUT_DIR}"/nonsyn_variants_MAF_*.tsv; do
  [ -e "$f" ] || { echo "No variants_MAF_*.tsv files found in ${EXTRACT_DIR}"; break; }

  # Make a safe tag for filenames (replace anything non-alnum with underscores)
  base="$(basename "$f" .tsv)"
  tag="$(printf '%s' "$base" | sed 's/[^A-Za-z0-9._-]/_/g')"   # e.g., variants_MAF_0.00_0.01

  echo ">>> Running PLINK --het for bin list: $f"
  plink \
    --bfile "${PLINK_CHR}" \
    --extract "$f" \
    --het \
    --allow-extra-chr \
    --out "${OUT_DIR}/het_${tag}"
done

SYN=outputs/allsyn

for f in "${OUT_DIR}"/syn_variants_MAF_*.tsv; do
  [ -e "$f" ] || { echo "No variants_MAF_*.tsv files found in ${EXTRACT_DIR}"; break; }

  # Make a safe tag for filenames (replace anything non-alnum with underscores)
  base="$(basename "$f" .tsv)"
  tag="$(printf '%s' "$base" | sed 's/[^A-Za-z0-9._-]/_/g')"   # e.g., variants_MAF_0.00_0.01

  echo ">>> Running PLINK --het for bin list: $f"
  plink \
    --bfile "${PLINK_CHR}" \
    --extract "$f" \
    --het \
    --allow-extra-chr \
    --out "${OUT_DIR}/het_${tag}"
done

