#!/bin/bash

input_file="ncbi_sequences.txt"
output_file="ncbi_sequences_upper_taxainfo.txt"

awk '
{
    # {} -> fro each line of the file:
    # Convert the sequence to uppercase
    $2 = toupper($2);
    # Extract the part of the line up to "taxon:tax_id" and remove everything else after -> $0 match the single line
    match($0, /(taxon:[^ ]+)/, taxon_match);
    # Keep the matched taxon code
    taxon_code = taxon_match[1];
    # Remove everything after the taxon code
    sub(/taxon:.*/, " " taxon_code);
    # Print the processed line
    print;
}' "$input_file" > "$output_file"

echo "Processing complete. Output written to $output_file."
