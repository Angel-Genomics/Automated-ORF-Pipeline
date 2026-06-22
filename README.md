# Automated DNA/RNA Sequence Analysis & BLAST Pipeline

## Overview
This repository contains a modular, 3-step bioinformatics pipeline written in Python. It is designed to automate the process of analyzing raw FASTA sequences downloaded from NCBI, extracting Open Reading Frames (ORFs), filtering the data for the most biologically relevant proteins, and querying them against the NCBI BLAST database.

## Pipeline Workflow & Files

The pipeline consists of three main scripts that must be executed sequentially:

### 1. `rna.py` (Extraction & Statistics)
* **Input:** `sequence.fasta` (Raw nucleotide sequences from NCBI)
* **Function:** Analyzes the sequence to instantly calculate nucleotide frequencies, total length, and GC/AU content in the terminal. It then identifies all possible ORFs across all 3 forward and 3 reverse reading frames.
* **Output:** Saves the extracted data into `all_orfs_results.csv` containing: `Gene_ID`, `Strand`, `Frame`, `Protein_Length`, and `Sequence`.

### 2. `filter_orf.py` (Data Transformation)
* **Input:** `all_orfs_results.csv`
* **Function:** Uses the `pandas` library to clean and filter the raw ORFs. It groups the data by `Gene_ID` and isolates the longest, most probable protein sequence for each gene.
* **Output:** Generates a refined dataset named `longest_orf_per_gene.csv`.

### 3. `blast.py` (Online Alignment)
* **Input:** `longest_orf_per_gene.csv`
* **Function:** Connects directly to the NCBI QBLAST API. It automates the BLAST search for the filtered protein sequences to find significant alignments in the `nr` database.
* **Output:** Produces the final analytical report `blast_results.csv`, containing critical biological metrics: `Gene_ID`, `Length`, `Blast Result`, `E-Value`, `Identity`, and `Protein Sequence`.

## How to Use
1. Download your desired nucleotide sequences from NCBI and save them as `sequence.fasta` in the root directory.
2. Run `rna.py` to extract all possible ORFs.
3. Run `filter_orf.py` to keep only the longest ORFs for each gene.
4. Run `blast.py` to fetch alignment results from the NCBI servers.

## Technologies Used
* Python 3
* Biopython (SeqIO, NCBIWWW, NCBIXML)
* Pandas (Data manipulation and CSV handling)
