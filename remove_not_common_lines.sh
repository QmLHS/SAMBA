#!/bin/bash

# Check if both files are provided in the bash command to execute the file
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <fasta_file> <taxonomy_file>"
    exit 1
fi

fasta_file=$1
tax_file=$2

# Check if files exist in the folder where the script is executed
if [ ! -f "$fasta_file" ] || [ ! -f "$tax_file" ]; then
    echo "Error: One or both input files do not exist"
    exit 1
fi

# Create output filenames adding '_filtered' to the same
fasta_output="${fasta_file%.fasta}_filtered.fasta"
tax_output="${tax_file%.tsv}_filtered.tsv"

temp_dir=$(mktemp -d)

grep '^>' "$fasta_file" | sed 's/^>//' | sort | uniq > "$temp_dir/fasta_ids.txt"
cut -f1 "$tax_file" | sort | uniq > "$temp_dir/tax_ids.txt"
# comm -12 keeps only the 3 flag that is the 'both' flag
comm -12 "$temp_dir/fasta_ids.txt" "$temp_dir/tax_ids.txt" > "$temp_dir/common_ids.txt"

# Filter FASTA file
# Initialize output file
> "$fasta_output"
# Flag to control writing
write_sequence=0
current_id=""

# VISUAL EXPLANATION
# Input FASTA:
# >ID1                  # Does this line start with '>'? If yes, is ID1 in common_ids? If yes, write_sequence=1, write this line
# ATGC                  # Does this line start with '>'? If no, write_sequence still 1, write this line
# GGTT                  # Does this line start with '>'? If no, write_sequence still 1, write this line
# >ID2                  # Does this line start with '>'? If yes, Is ID2 in common_ids? If no, write_sequence=0, do not write this line
# CCCC                  # Does this line start with '>'? If no, write_sequence=0, skip this line
# GGGG                  # Does this line start with '>'? If no, write_sequence still 0, skip this line
# >ID3                  # Does this line start with '>'? If yes, is ID3 in common_ids? If yes, write_sequence=1, write this line
# AATT                  # Does this line start with '>'? If no, write_sequence still 1, write this line
while IFS= read -r line; do
    if [[ $line =~ ^\> ]]; then
        # Extract ID from header line
        current_id=$(echo "$line" | sed 's/^>//')
        # Check if ID is in common_ids
        if grep -q "^${current_id}$" "$temp_dir/common_ids.txt"; then
            write_sequence=1
            echo "$line" >> "$fasta_output"
        else
            write_sequence=0
        fi
    elif [ $write_sequence -eq 1 ]; then
        echo "$line" >> "$fasta_output"
    fi
done < "$fasta_file"

# Filter taxonomy file
while IFS= read -r line; do
    id=$(echo "$line" | cut -f1)
    if grep -q "^${id}$" "$temp_dir/common_ids.txt"; then
        echo "$line" >> "$tax_output"
    fi
done < "$tax_file"

# Visual output
fasta_orig=$(grep -c '^>' "$fasta_file")
fasta_filt=$(grep -c '^>' "$fasta_output")
tax_orig=$(wc -l < "$tax_file")
tax_filt=$(wc -l < "$tax_output")

echo "Summary:"
echo "FASTA file: $fasta_orig entries reduced to $fasta_filt"
echo "Taxonomy file: $tax_orig entries reduced to $tax_filt"
echo "Filtered files saved as:"
echo "  FASTA: $fasta_output"
echo "  Taxonomy: $tax_output"
rm -rf "$temp_dir"
