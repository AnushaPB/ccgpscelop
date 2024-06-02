# The following is from: https://github.com/ccgproject/ccgpWorkflow/blob/fc4d5c89e5f8c55f428a7c867aeeb29714880cc3/workflow/modules/ccgp/common.smk#L8

import pandas as pd
import statistics
import gzip

def calc_roh(roh, fai, output):
    # Read the .fai file as a table using pandas, separating columns by tabs.
    dffai = pd.read_table(fai, sep='\t', header=None)
    # Filter the data to include only chromosomes greater than 1Mbp.
    dffai = dffai[dffai[1] > 10000000]  
    # Extract the chromosome identifiers into a list.
    chroms = dffai[0].values

    # Calculate the total length of the selected chromosomes.
    glength = dffai[1].sum()

    # Read the .roh file as a table, separating columns by tabs.
    dfroh = pd.read_table(roh, sep='\t', header=0)
    # Filter the ROH data to include segments greater than 500kbp.
    dfroh = dfroh[dfroh['[6]Length (bp)'] > 500000]  
    # Further filter to include only ROH data on the chromosomes selected above.
    dfroh = dfroh[dfroh['[3]Chromosome'].isin(chroms)]

    # Group the ROH data by sample, summing the lengths of ROH, and divide by the total genome length.
    dffroh = dfroh.groupby(['[2]Sample']).agg({'[6]Length (bp)': "sum"}).div(glength)
    # Reset the index of the grouped dataframe.
    dffroh = dffroh.reset_index()

    # Check if the grouped dataframe is empty (no data meets criteria).
    if dffroh.empty:
        # Print a message indicating no data was found.
        print("No individuals found with ROH greater than 500kbp.")
        # Create empty output files.
        with open(output, 'w'):
            pass
        with open(output.replace(".froh", "_top.froh"), "w"):
            pass
        # Return False indicating the operation was not successful.
        return False

    # Write the ROH data to a .froh file.
    dffroh.to_csv(output, sep='\t', index=False, header=None)
    # Print the dataframe for inspection.
    print(dffroh)

    # Check if the number of individuals is less than 8.
    if len(dffroh) < 8:
        # If true, use all the individuals in the final output.
        combined_df = dffroh
    else:
        # Select the two individuals with the largest ROH.
        df_lg = dffroh.nlargest(2, '[6]Length (bp)')
        # Randomly select 8 additional individuals.
        random_rows = dffroh.sample(n=8)
        # Combine the largest ROH individuals and random samples.
        combined_df = pd.concat([df_lg, random_rows], axis=0)

    # Write the combined data (either top 2 and 8 random, or all if less than 8) to a "_top.froh" file.
    combined_df.to_csv(output.replace(".froh", "_top.froh"), sep='\t', index=False, header=False, columns=[combined_df.columns[0]])
    # Return True indicating the operation was successful.
    return True