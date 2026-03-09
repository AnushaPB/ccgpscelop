calculate_ld_for_chr() {
  local BFILE="$1"        # PLINK prefix (no .bed/.bim/.fam)
  local CHR="$2"          # Chromosome number or name (e.g., 1, 2, GL123456)
  local WINDOW_SIZE=500000   # 500 kb
  local NSNPS=100             # 100 SNPs
  local TMPDIR="ld_chr${CHR}_tmp"
  local OUTFILE="ld_chr${CHR}_summary.txt"

  mkdir -p "$TMPDIR"
  > "$OUTFILE"

  # Get SNPs for this chromosome: CHR, SNPID, POS
  awk -v chr="$CHR" '$1 == chr {print $1, $2, $4}' "${BFILE}.bim" > "${TMPDIR}/chr_snps.txt"

  MIN_POS=$(awk '{print $3}' "${TMPDIR}/chr_snps.txt" | sort -n | head -1)
  MAX_POS=$(awk '{print $3}' "${TMPDIR}/chr_snps.txt" | sort -n | tail -1)

  for ((START=MIN_POS; START<=MAX_POS; START+=WINDOW_SIZE)); do
    END=$((START + WINDOW_SIZE - 1))

    awk -v s="$START" -v e="$END" '$3 >= s && $3 <= e {print $2}' "${TMPDIR}/chr_snps.txt" > "${TMPDIR}/snps_in_window.txt"

    if [ $(wc -l < "${TMPDIR}/snps_in_window.txt") -ge $NSNPS ]; then
      # Sample NSNPS SNPs
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
        done < <(tail -n +2 "${TMPDIR}/window_ld.ld")  # Skip header
      fi

      if [ "$N" -gt 0 ]; then
        MEAN_R2=$(echo "$R2_SUM / $N" | bc -l)
      else
        MEAN_R2="NA"
      fi

      echo -e "${CHR}\t${START}\t${END}\t${MEAN_R2}" >> "$OUTFILE"
    fi
  done

  echo "✅ Done: Chromosome ${CHR} LD summary in $OUTFILE"
}


PLINK=../../data/ccgp_data/58-Sceloporus_complete_coords_annotated
calculate_ld_for_chr $PLINK 5
mv ld_chr5_summary.txt outputs/