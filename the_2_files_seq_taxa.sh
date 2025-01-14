#!/bin/bash

input_file="ncbi_sequences_upper_taxainfo.txt"
fasta_file="ncbi_sequences.fasta"
taxon_file="ids_and_taxon.txt"

awk '
{
    id = $1;
    sequence = $2;
    print ">" id "\n" sequence;
}' "$input_file" > "$fasta_file"

echo "FASTA file created: $fasta_file"

awk '
{
    id = $1;
    match($0, /(taxon:[^ ]+)/, taxon_match);
    taxon_code = taxon_match[1];
    print id "\t" taxon_code;
}' "$input_file" > "$taxon_file"

echo "Taxon file created: $taxon_file"
