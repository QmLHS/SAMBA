#!/bin/bash

esearch -db nucleotide -email christian.persico@unimib.it \
  -query '(COI-5P [ALL] OR COI [ALL] OR CO1 [ALL] OR COX1 [ALL] OR Cytochrome Oxidase I [ALL] OR Cytochrome c Oxidase Subunit 1 [ALL] OR COX-1 [ALL] OR Cytochrome Oxidase Subunit I [ALL] OR COI gene [ALL] OR COI region [ALL] OR mtDNA COI [ALL] OR Mitochondrial COI [ALL] OR Mitochondrial Cytochrome Oxidase I [ALL] OR COI barcode [ALL] OR DNA barcoding region [ALL] OR COI fragment [ALL] OR COI-5P region [ALL])' | \
  efetch -format uid > accession_list.txt


# Function to process IDs in batches
# Function to process IDs in batches
process_batch() {
    local start=$1
    local batch_size=$2
    local total=$3
     
    # Process a batch of UIDs
    sed -n "${start},$(($start + $batch_size - 1))p" accession_list.txt | \
    efetch -db nucleotide -format gbc | \
    xtract -pattern INSDSeq -sep '\n' -tab '\t' \
        -element INSDSeq_primary-accession \
        -element INSDSeq_sequence \
        -block "INSDSeq_feature-table" \
        -subset INSDFeature \
        -if INSDFeature_key \
        -equals "source" \
        -subset INSDFeature_quals \
        -subset INSDQualifier \
        -if INSDQualifier_name \
        -equals "db_xref" \
        -element INSDQualifier_value >> ncbi_sequences.txt
     
    echo "Processed entries ${start} to $(($start + batch_size - 1)) of ${total}"
}

total_sequences=$(wc -l < accession_list.txt)
batch_size=1000

for ((i=1; i <= total_sequences; i+=batch_size)); do
    process_batch $i $batch_size $total_sequences
    sleep 5
done
