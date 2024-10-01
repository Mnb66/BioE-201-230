# Repo for Week 4-5 Assignment

## Command for question 1
seq = "KVRMFTSELDIMLSVNGPADQIKYFCRHWT"

print(len(seq))

print(len(seq)*3 + 3)

**Result:** 30 93

## Command for question 2 (Redirecting Prodigal's output to standard output)
prodigal -i genome.fna -o /dev/stdout | grep -c "CDS"

`-------------------------------------
PRODIGAL v2.6.3 [February, 2016]
Univ of Tenn / Oak Ridge National Lab
Doug Hyatt, Loren Hauser, et al.
-------------------------------------
Request:  Single Genome, Phase:  Training
Reading in the sequence(s) to train...1042519 bp seq created, 41.31 pct GC
Locating all potential starts and stops...37422 nodes
Looking for GC bias in different frames...frame bias scores: 2.60 0.20 0.20
Building initial set of genes to train from...done!
Creating coding model and scoring nodes...done!
Examining upstream regions and training starts...done!
-------------------------------------
Request:  Single Genome, Phase:  Gene Finding
Finding genes in sequence #1 (1042519 bp)...done!
897`

## Result for question 3
The biggest output is: "Genome with the highest number of genes: ./ncbi_dataset/ncbi_dataset/data/GCA_000006745.1/GCA_000006745.1_ASM674v1_genomic.fna (3594 genes)"
please refer to [question3.sh](https://github.com/Mnb66/BioE-201-230/blob/main/Genome%20Annotation/question3.sh) for code and [q3_gene_counts.txt](https://github.com/Mnb66/BioE-201-230/blob/main/Genome%20Annotation/q3_gene_counts.txt) for output result.

## Result for question 4 
please refer to [question4.sh](https://github.com/Mnb66/BioE-201-230/blob/main/Genome%20Annotation/question4.sh) for code and [gene_counts_difference.txt](https://github.com/Mnb66/BioE-201-230/blob/main/Genome%20Annotation/gene_counts_difference.txt) for output result.

## Command for question 5
find . -name "*.gff" -exec grep -oP 'Name=\K[^;]+' {} + | sort -u | head -n 5

### The output is like:
./prokka_output_GCA_000006745.1_ASM674v1_genomic/prokka_genes_GCA_000006745.1_ASM674v1_genomic.gff:aaeA
./prokka_output_GCA_000006745.1_ASM674v1_genomic/prokka_genes_GCA_000006745.1_ASM674v1_genomic.gff:aat
./prokka_output_GCA_000006745.1_ASM674v1_genomic/prokka_genes_GCA_000006745.1_ASM674v1_genomic.gff:abgT_1
./prokka_output_GCA_000006745.1_ASM674v1_genomic/prokka_genes_GCA_000006745.1_ASM674v1_genomic.gff:abgT_2
./prokka_output_GCA_000006745.1_ASM674v1_genomic/prokka_genes_GCA_000006745.1_ASM674v1_genomic.gff:accA
