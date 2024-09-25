# LLM：ChatGPT-4o
# 提示词：“编写一个 Python 脚本，使用起始密码子（‘ATG’）和终止密码子（‘TAA’、‘TAG’、‘TGA’）在基因组序列中查找基因。”

import sys

def read_fasta(file_path):
    with open(file_path, 'r') as f:
        sequence = ""
        for line in f:
            if not line.startswith('>'):
                sequence += line.strip()
    return sequence

def find_orfs(sequence):
    start_codon = "ATG"
    stop_codons = ["TAA", "TAG", "TGA"]
    orfs = []

    for frame in range(3):
        start = None
        for i in range(frame, len(sequence), 3):
            codon = sequence[i:i+3]
            if start is None and codon == start_codon:
                start = i
            elif start is not None and codon in stop_codons:
                orfs.append(sequence[start:i+3])
                start = None

    return orfs

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python gene_finder.py <fasta_file>")
        sys.exit(1)

    fasta_file = sys.argv[1]
    sequence = read_fasta(fasta_file)
    orfs = find_orfs(sequence)

    for orf in orfs:
        print(f"ORF found: {orf}")
