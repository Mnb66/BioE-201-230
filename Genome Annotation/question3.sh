#!/bin/bash

# Initialize variables
max_genes=0
max_genome=""

# Find all .fna files in the current directory and subdirectories
FILES=$(find . -name "*.fna")

# Loop over each genome file
for FILE in $FILES; do
    # Define output file name based on the input file name
    OUTPUT="${FILE%.fna}.gff"
    
    # Run prodigal on the genome file
    prodigal -i "$FILE" -o "$OUTPUT" -f gff > /dev/null 2>&1
    
    # Count the number of genes by counting the lines with "CDS" in the GFF output
    num_genes=$(grep -c -w "CDS" "$OUTPUT")
    
    # Display the genome file name and its gene count
    echo "$FILE has $num_genes genes."
    
    # Update the maximum gene count and corresponding genome file if current count is higher
    if [ "$num_genes" -gt "$max_genes" ]; then
        max_genes="$num_genes"
        max_genome="$FILE"
    fi
done

# Display the genome with the highest number of genes
echo "Genome with the highest number of genes: $max_genome ($max_genes genes)"

