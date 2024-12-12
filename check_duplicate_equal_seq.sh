#!/bin/bash

# Each line of duplicated_bold_records.txt is a couple -> (number, SequenceID)
for comp in $(awk '{print $2}' duplicated_bold_records.txt); do
  # Extract all sequences corresponding to the SequenceID in a vector
  sequences=($(grep -w $comp bold_allrawSeqs.fasta -A 1 | grep -v '^>' | grep -v '^--'))

  # Compare all pairs of sequences
  len=${#sequences[@]}
  for ((i = 0; i < len; i++)); do
    for ((j = i + 1; j < len; j++)); do
      str1=${sequences[i]}
      str2=${sequences[j]}
      if [ "$str1" == "$str2" ]; then
        echo "$comp: Sequence $i and Sequence $j are equal"
      else
        echo "$comp: Sequence $i and Sequence $j are not equal"
      fi
    done
  done
done
