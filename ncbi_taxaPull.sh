#!/bin/bash

process_batch() {
    local start=$1
    local batch_size=$2
    local total=$3
    
    # Process a batch of UIDs
    sed -n "${start},$(($start + $batch_size - 1))p" accession_list.txt | \
    xargs -n1 | \
      efetch -db taxonomy -format xml | \
      xtract -pattern Taxon -tab "," \
        -first TaxId ScientificName \
        -group Taxon \
        -KING "(-)" -PHYL "(-)" -CLSS "(-)" -ORDR "(-)" -FMLY "(-)" -GNUS "(-)" -SPCS "(-)" \
        -block "*/Taxon" -match "Rank:superkingdom" -KING ScientificName \
        -block "*/Taxon" -match "Rank:kingdom" -KING ScientificName \
        -block "*/Taxon" -match "Rank:phylum" -PHYL ScientificName \
        -block "*/Taxon" -match "Rank:class" -CLSS ScientificName \
        -block "*/Taxon" -match "Rank:order" -ORDR ScientificName \
        -block "*/Taxon" -match "Rank:family" -FMLY ScientificName \
        -block "*/Taxon" -match "Rank:genus" -GNUS ScientificName \
        -block "*/Taxon" -match "Rank:species" -SPCS ScientificName \
        -group Taxon -tab "," \
        -element "&KING" "&PHYL" "&CLSS" "&ORDR" "&FMLY" "&GNUS" "&SPCS" >> taxa_info.txt
    echo "Processed entries ${start} to $(($start + batch_size - 1)) of ${total}"
}

total_sequences=$(wc -l < accession_list.txt)
batch_size=1000

for ((i=1; i <= total_sequences; i+=batch_size)); do
    process_batch $i $batch_size $total_sequences
    sleep 5
done
