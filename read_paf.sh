PAF=map_GCA_029215755.1_aSpeHam1.0.hap2_genomic_to_jordan-uni4378-mb-hirise-65qvq__12-23-2023__final_assembly_relabelled.paf
head -n 5 $PAF 


awk '{sum += $11} END {print sum}' $PAF

# 803403665

# Filter to chr6 and sum
grep chr6 $PAF | awk '{sum += $10} END {print sum}'
# 19733450

# Filter to all overlap > 1 Mb and sum
grep -E '^[^#]' $PAF | awk '$10 > 1000000' | awk '{sum += $10} END {print sum}'
# 600329068