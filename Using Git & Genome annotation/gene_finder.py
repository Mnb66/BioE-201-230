# LLM：ChatGPT-4o
# Prompt: "Write a Python script to find a gene in a genome sequence using a start codon ('ATG') and a stop codon ('TAA', 'TAG ', 'TGA') 
# to find genes in the genome sequence."

import sys
from Bio.Seq import Seq

# Function to read a FASTA file and extract the sequence
def read_fasta(file_path):
    with open(file_path, 'r') as f:
        sequence = ""
        for line in f:
            if not line.startswith('>'):  # Skip header lines
                sequence += line.strip()
    return sequence

# Function to find Open Reading Frames (ORFs) in all six reading frames
def find_orfs_in_six_frames(sequence):
    # Generate all six reading frames
    frames = [sequence, sequence[1:], sequence[2:]]  # Forward frames
    reverse_complement = str(Seq(sequence).reverse_complement())
    frames += [reverse_complement, reverse_complement[1:], reverse_complement[2:]]  # Reverse frames

    orfs = []
    for frame in frames:
        orfs += find_orfs(frame)
    return orfs

# Function to find ORFs in a single sequence
def find_orfs(sequence):
    start_codon = "ATG"
    stop_codons = ["TAA", "TAG", "TGA"]
    orfs = []

    # Iterate through the sequence to find ORFs
    for frame in range(len(sequence)-2):
        start = None
        for i in range(frame, len(sequence), 3):
            codon = sequence[i:i+3]
            if start is None and codon == start_codon:
                start = i  # Found start codon
            elif start is not None and codon in stop_codons:
                orfs.append(sequence[start:i+3])  # Found stop codon, add ORF
                start = None  # Reset start position for next ORF

    return orfs

# Function to translate ORFs to protein sequences
def translate_orfs(orfs):
    protein_strings = set()
    for orf in orfs:
        protein = Seq(orf).translate()
        protein_strings.add(protein.split('*')[0])  # Remove stop codon and add to set
    return protein_strings

# Main execution block
if __name__ == "__main__":
    # Check if correct number of arguments is provided
    if len(sys.argv) != 2:
        print("Usage: python gene_finder.py <fasta_file>")
        sys.exit(1)
    fasta_file = sys.argv[1]
    sequence = read_fasta(fasta_file)  # Read the FASTA file
    orfs = find_orfs_in_six_frames(sequence)  # Find ORFs in all six frames
    # Translate ORFs to protein sequences and print them

    # Rosalind problem 72
    protein_strings = translate_orfs(orfs)
    for protein in protein_strings:
        print(protein)
