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


# Command line to geerate job lists
paste a2b_pairs.txt b2a_pairs.txt | \
  awk '{printf "python hungarian_scipy.py %s %s --emax 1000000 -o %s_hungarian.csv --summary %s_hungarian.summary\n",$1,$2,$1,$1}' \
  > jobs.txt

  # Use GNU parallel to run the job.txt jobs in parallel
  parallel -j20 < job.txt 
