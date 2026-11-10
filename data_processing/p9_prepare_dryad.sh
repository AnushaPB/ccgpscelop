
#!/bin/bash

# Run from scripts/data_processing/
# Copy files for Dryad archive

DATA=../data
DRYAD="$DATA/dryad"
CCGP="$DATA/ccgp_data"
GENOME="$DATA/genome"

mkdir -p "$DRYAD"

# VCF
cp "$CCGP/58-Sceloporus_complete_coords_annotated_final_samples.vcf.gz" "$DRYAD/"

# PLINK datasets (BED, BIM, FAM)
cp "$CCGP/58-Sceloporus_complete_coords_annotated_chr."{bed,bim,fam} "$DRYAD/"
cp "$CCGP/58-Sceloporus_annotated_pruned_0.6_chr."{bed,bim,fam} "$DRYAD/"
cp "$CCGP/58-Sceloporus_pruned_0.6_thinned_10kb_chr."{bed,bim,fam} "$DRYAD/"

# Sampling coordinates and callable sites
cp "$CCGP/58-Sceloporus.coords.txt" "$DRYAD/"
cp "$CCGP/58-Sceloporus_callable_sites.bed" "$DRYAD/"

# Reference genome and annotation
cp "$GENOME/rSceOcc1.20250428.p_ctg_chr_names.fna" "$DRYAD/"
cp "$GENOME/rSceOcc1.20250428.p_ctg_chr_names.fna.fai" "$DRYAD/"
cp "$GENOME/annotation/complete.genomic.gff" "$DRYAD/"

# GEA genes with nonsynonymous variants
cp "../analysis/gea/outputs/nonsynonymous_gea_genes.csv" "$DRYAD/"

# Copy README file
cp "README_dryad.md" "$DRYAD/README.md"

# SUBMIT TO DRYAD ------------------------------------------------

DRYAD_DOI="doi:10.5061/dryad.v15dv42df"
VCF="58-Sceloporus_complete_coords_annotated_final_samples.vcf.gz"

# Load API credentials stored outside the Git repository
source "$HOME/.config/dryad/credentials"

# Request a fresh token and extract it from the JSON response
DRYAD_TOKEN=$(curl --fail-with-body --silent --show-error \
  --data-urlencode "client_id=$DRYAD_CLIENT_ID" \
  --data-urlencode "client_secret=$DRYAD_CLIENT_SECRET" \
  --data-urlencode "grant_type=client_credentials" \
  "https://datadryad.org/oauth/token" |
  python3 -c 'import json, sys; print(json.load(sys.stdin)["access_token"])')

unset DRYAD_CLIENT_ID DRYAD_CLIENT_SECRET

# Encode the DOI for the API URL
DOI_ENCODED=$(python3 -c 'import sys; from urllib.parse import quote; print(quote(sys.argv[1], safe=""))' "$DRYAD_DOI")

# Verify authentication before starting uploads
curl --fail-with-body --silent --show-error \
  -H "Authorization: Bearer $DRYAD_TOKEN" \
  "https://datadryad.org/api/v2/test"
printf '\n'

# Collect all files, placing the largest VCF last
mapfile -d '' FILES < <(find "$DRYAD" -maxdepth 1 -type f ! -name "$VCF" -print0 | sort -z)
FILES+=("$DRYAD/$VCF")

# Upload files sequentially
for DRYAD_FILE in "${FILES[@]}"; do
  FILE_NAME=$(basename "$DRYAD_FILE")
  NAME_ENCODED=$(python3 -c 'import sys; from urllib.parse import quote; print(quote(sys.argv[1], safe=""))' "$FILE_NAME")

  # Set the content type based on the file extension
  case "$FILE_NAME" in
    *.vcf.gz) CONTENT_TYPE="application/gzip" ;;
    *.txt|*.bed|*.bim|*.fam|*.fai|*.gff|*.md) CONTENT_TYPE="text/plain" ;;
    *) CONTENT_TYPE="application/octet-stream" ;;
  esac

  echo "Uploading: $FILE_NAME"
  ls -lh "$DRYAD_FILE"

  # Upload each file and stop if the upload fails
  curl --fail-with-body --progress-bar --retry 1 \
    --upload-file "$DRYAD_FILE" \
    -H "Authorization: Bearer $DRYAD_TOKEN" \
    -H "Content-Type: $CONTENT_TYPE" \
    -H "Accept: application/json" \
    --write-out '\nHTTP status: %{http_code}\nElapsed: %{time_total} seconds\nUploaded: %{size_upload} bytes\n' \
    "https://datadryad.org/api/v2/datasets/$DOI_ENCODED/files/$NAME_ENCODED"

  printf '\nFinished: %s\n\n' "$FILE_NAME"
done

unset DRYAD_TOKEN
echo "All files uploaded successfully."