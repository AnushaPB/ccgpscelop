import pandas as pd

total_callable_sites = 0

with open("9-Laterallus_callable_sites.bed", "r") as file:
    for line in file:
        parts = line.strip().split()
        start = int(parts[1])
        end = int(parts[2])
        total_callable_sites += (end - start)

blra_callable = total_callable_sites

print("Total callable sites Laterallus:", total_callable_sites)

total_callable_sites = 0

with open("9-Rallus_callable_sites.bed", "r") as file:
    for line in file:
        parts = line.strip().split()
        start = int(parts[1])
        end = int(parts[2])
        total_callable_sites += (end - start)

vira_callable = total_callable_sites

# number of sites in dataset with only SNPs
# TODO: switch to occur within script
# calculated with this code:
# grep -vc "^#" 9-Rallus_snpsonly.recode.vcf #33,112,649
# grep -vc "^#" 9-Laterallus_snpsonly.recode.vcf #19,363,267
with open('QC/9-Laterallus_snpsonly_nsites.txt', 'r') as file:
    blra_snps = int(file.readline())

with open('QC/9-Rallus_snpsonly_nsites.txt', 'r') as file:
    vira_snps = int(file.readline())

# postfilter number of sites
# calculated with this code:
# zgrep -vc "^#" 9-Rallus_strictfilter.vcf.gz #591,607
# zgrep -vc "^#" 9-Laterallus_strictfilter.vcf.gz #224,005
with open('QC/9-Laterallus_strictfilter_nsites.txt', 'r') as file:
    blra_snps_filter = int(file.readline())

with open('QC/9-Rallus_strictfilter_nsites.txt', 'r') as file:
    vira_snps_filter = int(file.readline())

# calculate number of snps excluded by filter
vira_dif = vira_snps - vira_snps_filter
blra_dif = blra_snps - blra_snps_filter

# calculate the number of callable sites post filter (for heterozygosity calcs)
vira_callable_filter = vira_callable - vira_dif
blra_callable_filter = blra_callable - blra_dif

# Create DataFrame
df = pd.DataFrame({
    'species': ['Rallus limicola', 'Laterallus jamaicensis'],
    'callable_sites': [vira_callable, blra_callable],
    'total_SNPs': [vira_snps, blra_snps],
    'filtered_SNPs': [vira_snps_filter, blra_snps_filter],
    'SNPs_excluded_by_filter': [vira_dif, blra_dif],
    'callable_sites_post_filter': [vira_callable_filter, blra_callable_filter]
})

df.to_csv("callable_counts.csv", index = False)
print(df)
