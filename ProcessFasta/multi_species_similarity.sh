#!/bin/bash 
: <<'END_COMMENT'
Steps:

Prepare a list of species (fasta files) in a directory.
For each species, create a Diamond database.
Then for each pair (i, j) where i != j, run the two Diamond searches and compute the average.

Steps to get species similarity using diamond (https://github.com/bbuchfink/diamond)
# ============================================

1. Precompute databases for each species.
For each unordered pair (i, j) with i != j, do:
a. Run Diamond: i vs j -> get top hits for each i, average the identities -> avg_i_j
b. Run Diamond: j vs i -> get top hits for each j, average the identities -> avg_j_i
c. Overall similarity for pair (i,j) = (avg_i_j + avg_j_i) / 2
Store the result in a matrix.

END_COMMENT

#1. Precompute diamond databases

for fasta in *.fasta; do
diamond makedb --in $fasta -d ${fasta%.fasta}_db
done

# 2. Create species pairs
for i in *.fasta; do
  for j in *.fasta; do
    if [ "$i" < "$j" ]; then
      echo "${i%.fasta} ${j%.fasta}"
    fi
  done
done > pairs.txt
