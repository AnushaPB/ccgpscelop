from collections import defaultdict
import pandas as pd

# Read BED file
bed_file = "../../data/ccgp_data/58-Sceloporus_callable_sites.bed"
scaffold_lengths = defaultdict(int)

# Calculate scaffold lengths
with open(bed_file, 'r') as f:
    for line in f:
        scaffold, start, end = line.strip().split()
        scaffold_lengths[scaffold] += int(end) - int(start)

# Create a DataFrame
df = pd.DataFrame(list(scaffold_lengths.items()), columns=['Scaffold', 'Length'])

# Sort by length
df = df.sort_values(by='Length', ascending=False)

# Get the top 10 longest scaffolds
top_10 = df.head(10)

# Print the sorted DataFrame
print(top_10)

# Write the top 10 to a new file
top_10[['Scaffold']].to_csv("outputs/top_10_scaffolds.txt", sep="\t", index=False, header=False)

# Calculate what percentage of the genome is covered by the top 10 scaffolds
total_length = df['Length'].sum()
top_10_length = top_10['Length'].sum()
percentage_covered = (top_10_length / total_length) * 100
print(f"Percentage of the genome covered by the top 10 scaffolds: {percentage_covered:.2f}%")