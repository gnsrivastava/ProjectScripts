#!/bin/bash

# Read species (strip blanks/comments/leading-trailing spaces)
mapfile -t species < <(sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//' \
                           -e '/^$/d' -e '/^#/d' species.txt)

: > a2b_pairs.txt
: > b2a_pairs.txt

# Generate unordered pairs (i<j): A_vs_B and B_vs_A
for ((i=0; i<${#species[@]}-1; i++)); do
  for ((j=i+1; j<${#species[@]}; j++)); do
    A=${species[i]}
    B=${species[j]}
    printf '%s_vs_%s.tsv\n' "$A" "$B" >> a2b_pairs.txt
    printf '%s_vs_%s.tsv\n' "$B" "$A" >> b2a_pairs.txt
  done
done

echo "Wrote a2b_pairs.txt ($(wc -l < a2b_pairs.txt) lines)"
echo "Wrote b2a_pairs.txt ($(wc -l < b2a_pairs.txt) lines)"
