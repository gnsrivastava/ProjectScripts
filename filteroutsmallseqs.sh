# This is to remove any protein sequence with less than 50 amino acids. This is to ensure that the sequences are not the artifacts or peptides. 
# PMID: 20428524; 
# https://doi.org/10.1038/nbt1267
# change the path and .faa according to requirements

#!/bin/bash
mkdir -p filtered
INPUTFolder="/path to fasta directory"
find $INPUTFolder -type f -name "*.faa" \
| parallel -j 20 '
  base=$(basename {} .faa);
  python - <<EOF
from Bio import SeqIO
infile = "{}"
outfile = "filtered/${base}.faa"
with open(outfile, "w") as out:
    SeqIO.write(
        (r for r in SeqIO.parse(infile, "fasta") if len(r.seq) >= 50),
        out,
        "fasta"
    )
EOF
