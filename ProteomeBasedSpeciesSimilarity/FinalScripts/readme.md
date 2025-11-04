# Steps to DIAMOND BLASTP based assignment of proteins in BVBRC to STRING 

## Prerequisites
1. diamond
2. Biopython
3. Python v13
4. Protein Fasta (STRINGdb12)


**Step 1: Perform all the steps in ```ProcessBVBRC``` to generate input files for ```DIAMOND BLASTP```** 

**Step 2: For STRING bacteria strains create BlastpDB using diamond** 
- Run this script within the folder containing all the fasta files of each bacterial strain in STRINGdb v12
```
# Location of the fasta files StringSeqsOrg
ls *.fa | parallel --jobs 64 ' 
  file={} 
    base_name=$(basename "$file" .fa | cut -d"." -f1)
    location_to_diamond_installation/diamond makedb --in "$file" -d "$base_name" 
'
mv *.dmnd ../StringSeqsDiamond
```
Step 3: Run diamond on these non-redundant sequences against STRING (Script: 

Step 4: Run Hungarian Assignment (Script:
