import pandas as pd 
def transcribe_dna(dna_seq):
    clean_seq = dna_seq.replace('\n','')
    rna_seq = clean_seq.replace('T','U')
    return rna_seq
def analysis_rna(rna_seq):
    counts = {
        "A" : 0,
        "U" : 0,
        "G" : 0,
        "C" : 0
    }
    for noucleotide in rna_seq:
        if noucleotide in counts:
            counts[noucleotide] += 1
    print('\n---RNA ANALYSIS RESULTS---\n')
    print(f'Total length: {(len(rna_seq))}')
    for key, value in counts.items():
        print(f"Nucleotide {key} : {value}")
    print("\n----------------------------------\n")
    count_A = counts['A']
    count_U = counts['U']
    count_G = counts['G']
    count_C = counts['C']
    total_length = len(rna_seq)
    if total_length > 1:
        AU_content = ((count_A + count_U)/ total_length) * 100
        GC_content = ((count_G + count_C)/ total_length) * 100
        print(f'AU Content: {AU_content:.2f}%')
        print(f'GC Conetnt: {GC_content:.2f}%')
    else:
        print('ERROR: The file is Empty!!!')
        return
def get_reverse_complement(dna_seq):
    complement_map = {'A' : 'T' , 'T' : 'A' , 'C' : 'G' , 'G' : 'C'}
    complement_seq = ""
    for nucleotide in dna_seq:
        if nucleotide in complement_map:
            complement_seq += complement_map[nucleotide]
        else:
            complement_seq += nucleotide
    revers_comp_seq = complement_seq[::-1]
    return revers_comp_seq
codon_table = {
        'AUA': 'I', 'AUC': 'I', 'AUU': 'I', 'AUG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACU': 'T',
        'AAC': 'N', 'AAU': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGU': 'S', 'AGA': 'R', 'AGG': 'R',
        'CUA': 'L', 'CUC': 'L', 'CUG': 'L', 'CUU': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCU': 'P',
        'CAC': 'H', 'CAU': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGU': 'R',
        'GUA': 'V', 'GUC': 'V', 'GUG': 'V', 'GUU': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCU': 'A',
        'GAC': 'D', 'GAU': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGU': 'G',
        'UCA': 'S', 'UCC': 'S', 'UCG': 'S', 'UCU': 'S',
        'UUC': 'F', 'UUU': 'F', 'UUA': 'L', 'UUG': 'L',
        'UAC': 'Y', 'UAU': 'Y', 'UAA': '*', 'UAG': '*',
        'UGC': 'C', 'UGU': 'C', 'UGA': '*', 'UGG': 'W',
    }
def find_all_orfs(dna_seq, gene_id, min_port_len=100):
    all_orfs = []
    rev_com_dna = get_reverse_complement(dna_seq)
    forwad_rna = transcribe_dna(dna_seq)
    revers_rna = transcribe_dna(rev_com_dna)
    strands = {'+' : forwad_rna, '-' : revers_rna}
    for strand_symbol, rna_seq in strands.items():
        for frame in range(3):
            is_coding = False
            current_protein = ""
            for i in range(frame, len(rna_seq) - 2, 3):
                codon = rna_seq[i:i+3]
                if codon == "AUG" and not is_coding:
                    is_coding = True
                    current_protein = "M"
                elif codon in ["UAA", "UGA", "UAG"] and is_coding:
                    is_coding = False
                    if len(current_protein) >= min_port_len:
                        all_orfs.append({
                            'Gene_ID': gene_id,
                            'Strand': strand_symbol,
                            'Frame': frame,
                            'Protein_Length': len(current_protein),
                            'Sequence': current_protein
                        })
                    current_protein = ""
                elif is_coding:
                    amino_acid = codon_table.get(codon, "X")
                    current_protein += amino_acid
    return all_orfs
if __name__ == '__main__':
    file_name = "sequence.fasta"
    current_name = ""
    sequences = {}
    all_results = []
    with open(file_name, "r") as file:
        for line in file:
            clean_line = line.strip()
            if len(clean_line) == 0:
                continue
            if clean_line.startswith(">"):
                current_name = clean_line[1:]
                sequences[current_name] = ""
            else:
                sequences[current_name] += clean_line
    print(f"---Successfully Read The {file_name}---\n")
    print(f"There is {len(sequences)} Sequences Found in The {file_name}\n")
    for gene_name, dna_seq in sequences.items():
        gene_id = gene_name.split()[0]
        print('='*50)
        print(f'Processing Gene Name: {gene_id}')
        print(f'Gene Full Description:\n {gene_name}')
        current_rna = transcribe_dna(dna_seq)
        analysis_rna(current_rna)
        gene_orfs = find_all_orfs(dna_seq, gene_id)
        all_results.extend(gene_orfs)
        print(f'Found {len(gene_orfs)} ORFs For {gene_id}')
    print('\n' + '='*50)
    print('---Saving Results To CSV---')
    df = pd.DataFrame(all_results)
    output_csv = 'all_orfs_results.csv'
    df.to_csv(output_csv, index=False)
    print(f'---Seccessfully Saved All {len(df)} ORFs To {output_csv}---') 