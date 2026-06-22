import pandas as pd
input_csv_file = 'all_orfs_results.CSV'
output_csv_file = 'longest_orf_per_gene.csv'
try:
    print(f'Reading Data From {input_csv_file} ...')
    df = pd.read_csv(input_csv_file)
    if df.empty:
        print(f"The {input_csv_file} is EMPTY!!!")
    else:
        longest_orf_index = df.groupby('Gene_ID')['Protein_Length'].idxmax()
        longest_orf_df_per_gene = df.loc[longest_orf_index]
        longest_orf_df_per_gene.to_csv(output_csv_file, index=False)
        print(f"---\nSuccessfully Found Longest ORF in the {input_csv_file} And Save it To {output_csv_file}.")
except FileNotFoundError:
    print(f"ERROR: The File {input_csv_file} Was Not Found. Please Make Sure it is in The Same Directory!!!")
except KeyError:
    print("ERROR: The CSV File Must Contain a Column Named 'Gene_ID' And 'Protein_Length'. Please Check The File.")
except Exception as e:
    print(f"An Unexpected ERROR Occurred: {e}")