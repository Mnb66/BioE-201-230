# LLM：ChatGPT-4o
# Prompt: "Implement numpy to speed up the gene finder."

import sys
from Bio.Seq import Seq
from tqdm import tqdm
import numpy as np

# Function to read a FASTA file and extract the sequence
def read_fasta(file_path):
    """
    Reads a FASTA file and extracts the sequence.
    
    Args:
    file_path (str): Path to the FASTA file.
    
    Returns:
    str: The extracted DNA sequence.
    """
    with open(file_path, 'r') as f:
        sequence = ""
        for line in f:
            if not line.startswith('>'):  # Skip header lines
                sequence += line.strip()
    return sequence

# Function to find Open Reading Frames (ORFs) in all six reading frames using NumPy
def find_orfs_in_six_frames(sequence, min_length=100):
    """
    Finds ORFs in all six reading frames of a DNA sequence.
    
    Args:
    sequence (str): The input DNA sequence.
    
    Returns:
    list: A list of all ORFs found in the six frames.
    """
    orfs = []
    frames = [0, 1, 2]
    # Forward frames
    for frame in frames:
        orfs += find_orfs_numpy(sequence, frame, min_length)
    # Reverse frames
    reverse_complement = str(Seq(sequence).reverse_complement())
    for frame in frames:
        orfs += find_orfs_numpy(reverse_complement, frame, min_length)
    return orfs

# Function to check if a sequence has a Shine-Dalgarno sequence
def has_shine_dalgarno(upstream_seq):
    sd_sequence = "AGGAGG"
    return sd_sequence in upstream_seq

# Function to find ORFs in a single sequence using NumPy
def find_orfs_numpy(sequence, frame=0, min_length=100):
    """
    Finds ORFs in a single frame of a DNA sequence using NumPy for efficiency.
    
    Args:
    sequence (str): The input DNA sequence.
    frame (int): The reading frame (0, 1, or 2).
    
    Returns:
    list: A list of ORFs found in the specified frame.
    """
    stop_codons = ['TAA', 'TAG', 'TGA']
    sequence = sequence[frame:]
    chars = np.array(list(sequence))
    n_codons = len(chars) // 3
    if n_codons == 0:
        return []
    codons = chars[:n_codons*3].reshape(-1, 3)
    codon_strings = np.array([''.join(codon) for codon in codons])
    start_indices = np.where(codon_strings == 'ATG')[0]
    stop_indices = np.where(np.isin(codon_strings, stop_codons))[0]

    orfs = []
    if len(start_indices) == 0 or len(stop_indices) == 0:
        return orfs

    # For each start_index, find the index in stop_indices where stop_idx >= start_idx
    idx_positions = np.searchsorted(stop_indices, start_indices, side='left')

    for i, start_idx in enumerate(start_indices):
        stop_idx_pos = idx_positions[i]
        # Ensure that the stop index is within bounds
        if stop_idx_pos < len(stop_indices):
            stop_idx = stop_indices[stop_idx_pos]
            if stop_idx >= start_idx:
                orf_length = stop_idx - start_idx + 1
                if orf_length >= min_length: # ORF min_length
                    # Check if the ORF has a Shine-Dalgarno sequence
                    upstream_start = max(0, start_idx * 3 - 20)
                    upstream_seq = sequence[upstream_start:start_idx * 3]
                    if has_shine_dalgarno(upstream_seq):
                        orf_codons = codon_strings[start_idx:stop_idx+1]
                        orf_seq = ''.join(orf_codons)
                        orfs.append(orf_seq)
    return orfs

# Function to translate ORFs to protein sequences
def translate_orfs(orfs):
    """
    Translates a list of ORFs to protein sequences.
    
    Args:
    orfs (list): A list of DNA sequences representing ORFs.
    
    Returns:
    set: A set of unique protein sequences.
    """
    protein_strings = set()
    for orf in tqdm(orfs, desc="Translating ORFs"):
        protein = Seq(orf).translate(to_stop=True)
        protein_strings.add(str(protein))
    return protein_strings

if __name__ == "__main__":
    # Check if the correct number of command-line arguments is provided
    if len(sys.argv) != 2:
        print("Usage: python gene_finder.py <fasta_file>")
        sys.exit(1)
    
    # Get the FASTA file path from command-line argument
    fasta_file = sys.argv[1]
    
    # Read the FASTA file and extract the sequence
    sequence = read_fasta(fasta_file)
    
    # Find ORFs in all six frames
    orfs = find_orfs_in_six_frames(sequence)
    
    # Translate ORFs to protein sequences
    protein_strings = translate_orfs(orfs)
    
    # Print each unique protein sequence
    # for protein in protein_strings:
    #     print(protein)
    
    print('Job for {} finished'.format(fasta_file))
