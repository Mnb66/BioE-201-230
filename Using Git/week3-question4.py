import sys
from Bio.Seq import Seq
from tqdm import tqdm
import numpy as np

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

def find_orfs_numpy_q_4(sequence,frame=0):
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
                orf_codons = codon_strings[start_idx:stop_idx+1]
                orf_seq = ''.join(orf_codons)
                orfs.append(orf_seq)
    return orfs

def find_orfs_in_six_frames(sequence, func, min_length=None):
    orfs = []
    frames = [0, 1, 2]
    for frame in frames:
        if min_length is not None:
            orfs += func(sequence, frame, min_length)
        else:
            orfs += func(sequence, frame)
    reverse_complement = str(Seq(sequence).reverse_complement())
    for frame in frames:
        if min_length is not None:
            orfs += func(reverse_complement, frame, min_length)
        else:
            orfs += func(reverse_complement, frame)
    return orfs

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python week3-question4.py <fasta_file>")
        sys.exit(1)
    fasta_file = sys.argv[1]
    sequence = read_fasta(fasta_file)
    orfs = find_orfs_in_six_frames(sequence, find_orfs_numpy_q_4)
    protein_strings = translate_orfs(orfs)
    print('Job for file {} finished'.format(fasta_file))
    # for protein in protein_strings:
    #     print(protein)
