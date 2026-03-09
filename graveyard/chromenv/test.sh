calculate_ld_for_chr() {
  local BFILE="$1"
  local CHR="$2"
  local WINDOW_SIZE=500000
  local NSNPS=100
  local TMPDIR="ld_chr${CHR}_tmp"
  local OUTFILE="ld_chr${CHR}_summary.txt"
  local SNP_DENSITY_FILE="ld_chr${CHR}_snp_density.txt"

  mkdir -p "$TMPDIR"
  > "$OUTFILE"
  > "$SNP_DENSITY_FILE"

  # 1. Extract SNPs (already good)
  awk -v chr="$CHR" '$1 == chr {print $1, $2, $4}' "${BFILE}.bim" > "${TMPDIR}/chr_snps.txt"

  MIN_POS=$(awk '{print $3}' "${TMPDIR}/chr_snps.txt" | sort -n | head -1)
  MAX_POS=$(awk '{print $3}' "${TMPDIR}/chr_snps.txt" | sort -n | tail -1)

  # 2. FIRST LOOP: LD (r²) calculation
  for ((START=MIN_POS; START<=MAX_POS; START+=WINDOW_SIZE)); do
    END=$((START + WINDOW_SIZE - 1))

    awk -v s="$START" -v e="$END" '$3 >= s && $3 <= e {print $2}' "${TMPDIR}/chr_snps.txt" > "${TMPDIR}/snps_in_window.txt"

    if [ $(wc -l < "${TMPDIR}/snps_in_window.txt") -ge $NSNPS ]; then
      # Sample SNPs
      shuf -n $NSNPS "${TMPDIR}/snps_in_window.txt" > "${TMPDIR}/sampled_snps.txt"

      # Calculate pairwise LD among sampled SNPs
      plink --bfile "$BFILE" \
            --extract "${TMPDIR}/sampled_snps.txt" \
            --r2 --ld-window 99999 --ld-window-kb 50000 --ld-window-r2 0 \
            --allow-extra-chr \
            --out "${TMPDIR}/window_ld" \
            --silent > /dev/null 2>&1

      R2_SUM=0
      N=0

      if [ -f "${TMPDIR}/window_ld.ld" ]; then
        while read -r line; do
          R2=$(echo "$line" | awk '{print $7}')
          if [[ $R2 =~ ^[0-9.]+$ ]]; then
            R2_SUM=$(echo "$R2_SUM + $R2" | bc -l)
            N=$((N + 1))
          fi
        done < <(tail -n +2 "${TMPDIR}/window_ld.ld")
      fi

      if [ "$N" -gt 0 ]; then
        MEAN_R2=$(echo "$R2_SUM / $N" | bc -l)
      else
        MEAN_R2="NA"
      fi

      echo -e "${CHR}\t${START}\t${END}\t${MEAN_R2}" >> "$OUTFILE"
    fi
  done

  # 3. SECOND LOOP: SNP density (after LD done)
  for ((START=MIN_POS; START<=MAX_POS; START+=WINDOW_SIZE)); do
    END=$((START + WINDOW_SIZE - 1))

    # Count SNPs in this window
    SNP_COUNT=$(awk -v s="$START" -v e="$END" '$3 >= s && $3 <= e' "${TMPDIR}/chr_snps.txt" | wc -l)

    echo -e "${CHR}\t${START}\t${END}\t${SNP_COUNT}" >> "$SNP_DENSITY_FILE"
  done

  echo "✅ Done: Chromosome ${CHR} LD and SNP density summary in $OUTFILE and $SNP_DENSITY_FILE"
}

PLINK=../../data/ccgp_data/58-Sceloporus_complete_coords_annotated
calculate_ld_for_chr $PLINK 5
mv ld_chr5_summary.txt outputs/
mv ld_chr5_snp_density.txt outputs/
