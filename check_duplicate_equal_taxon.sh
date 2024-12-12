#!/bin/bash

# Each line of duplicated_bold_records.txt is a couple -> (number, SequenceID)
# read line by line the file given;
# IFS= ensure no leading/trailing whitespace is removed from the lines because this way we are telling bash that there is no field separator;
# -r option do not let backslashes being interpreted as escape char.
while IFS= read -r line; do

  # Store second field: SequenceID
  comp=$(echo "$line" | awk '{print $2}')

  # Extract all the headers in an array, handling white spaces
  # mapfile: read output and put it inside array headers; -t prevent newlines from inclusion;
  # grep -w:  search for whole words |
  #    while: loop for each line in file read and for each header |
  #    cut: take the taxa reference being the second and so on element divided by ; |
  #    sed: remove leading and trainling white spaces;
  mapfile -t headers < <(grep -w "$comp" bold_allrawSeqs.fasta | while read -r header; do
    echo "$header" | cut -d ';' -f 2- | sed 's/^ *//; s/ *$//'
  done)

  # Compare all pairs of sequences
  len=${#headers[@]}
  for ((i = 0; i < len; i++)); do
    for ((j = i + 1; j < len; j++)); do
      taxon1="${headers[i]}"
      taxon2="${headers[j]}"
      if [ "$taxon1" == "$taxon2" ]; then
        echo "$comp: Taxa $i and Taxa $j are equal"
      else
        echo "Taxa $i: $taxon1"
        echo "Taxa $j: $taxon2"
        echo "$comp: Taxa $i and Taxa $j are not equal"
      fi
    done
  done
done < duplicated_bold_records.txt
