# !/bin/bash

# Find all .fna files in the current directory and subdirectories
FILES=$(find . -name "*.fna")

# Loop over each genome file
for FILE in $FILES; do
    # Define output file name based on the input file name
    OUTPUT_Prodgal="${FILE%.fna}_prodgal.gff"

    OUTPUT_Prokka="prokka_output_$(basename "$FILE" .fna)"
    prefix_Prokka="prokka_genes_$(basename "$FILE" .fna)"
    # Run prodigal on the genome file
    prodigal -i "$FILE" -o "$OUTPUT_Prodgal" -f gff > /dev/null 2>&1

    prokka --outdir "$OUTPUT_Prokka" --prefix "$prefix_Prokka" --quiet "$FILE" > /dev/null 2>&1

    num_genes_prokka=$(grep -c "CDS" "$OUTPUT_Prokka/$prefix_Prokka.gbk")
    
    # Count the number of genes by counting the lines with "CDS" in the GFF output
    num_genes_prodigal=$(grep -c -w "CDS" "$OUTPUT_Prodgal")
    
    # Display the genome file name and its gene count
    echo "$FILE has $num_genes_prodigal genes using Prodigal." >> gene_counts_prodigal.txt
    echo "$FILE has $num_genes_prokka genes using Prokka." >> gene_counts_prokka.txt

    # show the difference between the two
    echo "The gene num difference between using Prodigal and Prokka in $FILE is $((num_genes_prokka - num_genes_prodigal))" >> gene_counts_difference.txt
done