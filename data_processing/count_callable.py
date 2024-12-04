import pandas as pd

total_callable_sites = 0

with open("../data/ccgp_data/58-Sceloporus_callable_sites.bed", "r") as file:
    for line in file:
        parts = line.strip().split()
        start = int(parts[1])
        end = int(parts[2])
        total_callable_sites += (end - start)

callable = total_callable_sites

# Create DataFrame
df = pd.DataFrame({
    'callable_sites': [callable]
})

df.to_csv("callable_counts.csv", index = False)
print(df)
