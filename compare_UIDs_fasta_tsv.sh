#!/bin/bash

# Check if both files are provided in the input: first the FASTA then the tsv
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <fasta_file> <taxonomy_file>"
    exit 1
fi

fasta_file=$1
tax_file=$2

# Check if files exist in the directory where the script is executed
if [ ! -f "$fasta_file" ] || [ ! -f "$tax_file" ]; then
    echo "Error: One or both input files do not exist"
    exit 1
fi

fasta_ids=$(grep '^>' "$fasta_file" | sed 's/^>//' | sort | uniq)
tax_ids=$(cut -f1 "$tax_file" | sort | uniq)

temp_dir=$(mktemp -d)
echo "$fasta_ids" > "$temp_dir/fasta_ids.txt"
echo "$tax_ids" > "$temp_dir/tax_ids.txt"

fasta_only=$(comm -23 "$temp_dir/fasta_ids.txt" "$temp_dir/tax_ids.txt")
tax_only=$(comm -13 "$temp_dir/fasta_ids.txt" "$temp_dir/tax_ids.txt")

fasta_only_count=$(echo "$fasta_only" | grep -c "^" || true)
tax_only_count=$(echo "$tax_only" | grep -c "^" || true)

echo "IDs in FASTA file but not in taxonomy file: $fasta_only_count"
echo "IDs in taxonomy file but not in FASTA file: $tax_only_count"

rm -rf "$temp_dir"
