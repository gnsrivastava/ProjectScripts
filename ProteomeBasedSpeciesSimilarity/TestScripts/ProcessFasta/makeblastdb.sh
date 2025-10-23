ls *.faa | parallel --jobs 20 '
file={}
  base_name=$(basename "$file" .faa | cut -d"_" -f1)
  mkdir -p "$base_name"
  $makeblastdb -in "$file" -dbtype prot -out "$base_name/${base_name}_db"
'

# Script for Diamond makedb
mkdir OUTDIR
ls INPUTDIR/*.faa | parallel --jobs 20 '
  file={}
  base_name=$(basename "$file" .faa | cut -d"_" -f1)
  /INSTALLATIONLOCATION/diamond makedb --in {} -d ./OUTDIR/$base_name
'
