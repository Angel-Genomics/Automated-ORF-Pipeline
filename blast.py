import pandas as pd
from Bio.Blast import NCBIWWW, NCBIXML
import time
import os
input_file = "longest_orf_per_gene.csv"
output_file = "blast_results.csv"
df_input = pd.read_csv(input_file)
total_gene = len(df_input)
processed_genes = []
if os.path.exists(output_file):
    df_output = pd.read_csv(output_file)
    if "Gene_ID" in df_output.columns:
        processed_genes = df_output["Gene_ID"].tolist()
        print(f"Found {len(processed_genes)} Previously Precessed Genes.Resuming...")
write_header = not os.path.exists(output_file)
for index, row in df_input.iterrows():
    gene_id = row["Gene_ID"]
    if gene_id in processed_genes:
        print(f"Skipping {gene_id} ({index + 1} of {total_gene}). Already processed.")
        continue
    protein_sequence = str(row['Sequence']).strip()
    protein_length = str(row['Protein_Length']).strip()
    print("="*50)
    print(f"[*] Processing Gene: {gene_id} ({index + 1} of {total_gene})")
    print("Sending sequence to NCBI BLAST. This may take a few minutes.")
    try:
        blast_results = NCBIWWW.qblast("blastp", "nr", protein_sequence)
        blast_record = NCBIXML.read(blast_results)
        best_alignment = "No Significant Similarity Found!!!"
        e_value = "None"
        identity_str = 'None'
        if blast_record.alignments:
            alignment = blast_record.alignments[0]
            hsp = alignment.hsps[0]
            best_alignment = alignment.title
            e_value = hsp.expect
            identity_percent = (hsp.identities / hsp.align_length) * 100
            identity_str = f"{identity_percent:.2f}%"
            print(f"Results: {best_alignment[:60]}... | E_vlaue: {e_value}")
        new_results = pd.DataFrame([{
            'Gene_ID': gene_id,
            'Length': protein_length,
            'Blast Result': best_alignment,
            'E-Value': e_value,
            'Identity': identity_str,
            'Protein Sequence': protein_sequence
        }])
        new_results.to_csv(output_file, mode="a", header=write_header, index=False)
        write_header = False
        time.sleep(10)
    except Exception as e:
        print(f"An Unexpected ERROR for {gene_id} eccurred: {e}")
        print("Network ERROR. Stopping the script. Run again to resume from this Gene.")
        break