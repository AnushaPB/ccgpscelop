calculate_ld_for_chr() {
  local BFILE="$1"        # PLINK prefix (no .bed/.bim/.fam)
  local CHR="$2"          # Chromosome number or name (e.g., 1, 2, GL123456)
  local WINDOW_SIZE=100000
  local NSNPS=10
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
      shuf -n $NSNPS "${TMPDIR}/snps_in_window.txt" > "${TMPDIR}/sampled_snps.txt"
      R2_SUM=0
      N=0

      mapfile -t SNPS < "${TMPDIR}/sampled_snps.txt"
      for ((i = 0; i < ${#SNPS[@]}; i++)); do
        for ((j = i + 1; j < ${#SNPS[@]}; j++)); do
          SNP1=${SNPS[i]}
          SNP2=${SNPS[j]}
          plink --bfile "$BFILE" --ld "$SNP1" "$SNP2" \
                --allow-extra-chr --out "${TMPDIR}/pair_ld" \
                --noweb --silent > /dev/null 2>&1
          R2=$(grep "R-sq =" "${TMPDIR}/pair_ld.log" | awk '{print $4}')
          if [[ $R2 =~ ^[0-9.]+$ ]]; then
            R2_SUM=$(echo "$R2_SUM + $R2" | bc -l)
            N=$((N + 1))
          fi
        done
      done

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