import pandas as pd

bed_file = "../../data/ccgp_data/58-Sceloporus_callable_sites.bed"
bed_df = pd.read_csv(bed_file, sep="\t", header=None, names=["scaffold", "start", "end"])

# Calculate lengths for each row, sum by scaffold, and reset the index
df = (bed_df.assign(length=bed_df["end"] - bed_df["start"])
      .groupby("scaffold")["length"].sum()
      .reset_index())

# Sort by length
df = df.sort_values(by="length", ascending=False)

# Get the top 10 longest scaffolds
top_10 = df.head(10)
print(top_10)

# Write the top 10 to a file which will be read in for SMC++
top_10[['scaffold']].to_csv("outputs/top_10_scaffolds.txt", sep="\t", index=False, header=False)

# Calculate what percentage of the genome is covered by the top 10 scaffolds
total_length = df['length'].sum()
top_10_length = top_10['length'].sum()
percentage_covered = (top_10_length / total_length) * 100
print(f"Percentage of the genome covered by the top 10 scaffolds: {percentage_covered:.2f}%")