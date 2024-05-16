import pandas as pd

total_callable_sites = 0

with open("../raw_data/58-Sceloporus_callable_sites.bed", "r") as file:
    for line in file:
        parts = line.strip().split()
        start = int(parts[1])
        end = int(parts[2])
        total_callable_sites += (end - start)

callable = total_callable_sites

total_callable_sites = 0

# number of sites in dataset with only SNPs
# TODO: switch to occur within script
# calculated with this code:
# grep -vc "^#" $RAW_DATA > 58-Sceloporus_snpsonly_nsites.txt
with open('58-Sceloporus_allsnps_nsites.txt', 'r') as file:
    snps = int(file.readline())

# postfilter number of sites
# calculated with this code:
# zgrep -vc "^#" 58-Sceloporus_maf05_minDP5_maxDP50_rmsamp60_mm80_rmsamp20.vcf.gz > 58-Sceloporus_postfiltersnps_nsites.txt
with open('58-Sceloporus_postfiltersnps_nsites.txt', 'r') as file:
    snps_filter = int(file.readline())

# calculate number of snps excluded by filter
dif = snps - snps_filter

# calculate the number of callable sites post filter (for heterozygosity calcs)
callable_filter = callable - dif

# Create DataFrame
df = pd.DataFrame({
    'callable_sites': [callable],
    'total_SNPs': [snps],
    'filtered_SNPs': [snps_filter],
    'SNPs_excluded_by_filter': [dif],
    'callable_sites_post_filter': [callable_filter]
})

df.to_csv("callable_counts.csv", index = False)
print(df)
