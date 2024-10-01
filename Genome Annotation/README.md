# Repo for Week 4-5 Assignment

## Command for question 1
seq = "KVRMFTSELDIMLSVNGPADQIKYFCRHWT"
print(len(seq))
print(len(seq)*3 + 3)
**Result:** 30 93

## Command for question 2 (Redirecting Prodigal's output to standard output)
prodigal -i genome.fna -o /dev/stdout | grep -c "CDS"

## Command for question 5
find . -name "*.gff" -exec grep -oP 'Name=\K[^;]+' {} + | sort -u | head -n 5

### The output is like:
./prokka_output_GCA_000006745.1_ASM674v1_genomic/prokka_genes_GCA_000006745.1_ASM674v1_genomic.gff:aaeA
./prokka_output_GCA_000006745.1_ASM674v1_genomic/prokka_genes_GCA_000006745.1_ASM674v1_genomic.gff:aat
./prokka_output_GCA_000006745.1_ASM674v1_genomic/prokka_genes_GCA_000006745.1_ASM674v1_genomic.gff:abgT_1
./prokka_output_GCA_000006745.1_ASM674v1_genomic/prokka_genes_GCA_000006745.1_ASM674v1_genomic.gff:abgT_2
./prokka_output_GCA_000006745.1_ASM674v1_genomic/prokka_genes_GCA_000006745.1_ASM674v1_genomic.gff:accA
