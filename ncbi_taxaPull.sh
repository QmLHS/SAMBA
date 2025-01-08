#!/bin/bash

cat accession_list.txt | \
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
