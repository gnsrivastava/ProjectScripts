: <<'END_COMMENT'
The script takes in the list of fasta files in original bvbrc and the splits the lists into chunks of 2000 files per chunk
for each chunk then a GNU parallel based SLRUM job submisison script is created in "slurm_scripts" 
The scripts in the slurm_scripts can then be submitted as follows

for file in slurm_scripts/*.sh
  do
  sbatch $file;
  sleep 5;
done;

END_COMMENT

#!/bin/bash

# === Configuration ===
IDS="lst"              # master input file (list of .faa files)
CHUNK_DIR="./chunks"       # where chunk files will live
SCRIPT_DIR="./slurm_scripts" # where per-chunk Slurm scripts will live
LINES_PER_CHUNK=2000        # **Adjust this** to control number of files per chunk
JOBS_PER_CHUNK=64          # **Adjust this** number of concurrent tasks per chunk (max 'n' in #SBATCH -n)

SLRUM_RUN_DIR=$( pwd )
INPUT_DIR="fastaBVBRC_org"
OUTPUT_DIR="filtered"

mkdir -p "$CHUNK_DIR" "$SCRIPT_DIR" "$OUTPUT_DIR"

# Clean old files + scripts
rm -f "$IDS" "$CHUNK_DIR"/part_* "$SCRIPT_DIR"/chunk_*.sh

# 1. === Generate the list of FASTA files to process ===
# Get relative paths of all .faa files and write them to the IDS file
find ./"$INPUT_DIR" -type f -name "*.faa" > "$IDS"

if [ ! -s "$IDS" ]; then
    echo "Error: No *.faa files found in $INPUT_DIR. Exiting."
    exit 1
fi

# 2. === Split into chunks ===
split -l "$LINES_PER_CHUNK" --numeric-suffixes=0 --suffix-length=5 \
  "$IDS" "$CHUNK_DIR/part_"

# 3. === Make a Slurm script for each chunk ===
for chunk in "$CHUNK_DIR"/part_*; do
  base=$(basename "$chunk")
  jobname="faa_filter_${base}"

  cat > "$SCRIPT_DIR/${base}.sh" <<EOF
#!/bin/bash

#SBATCH -J $jobname
#SBATCH -o $jobname.%j.out
#SBATCH -e $jobname.%j.err
#SBATCH -p workq
#SBATCH -t 00:10:00  # Adjust time as needed
#SBATCH -N 1
#SBATCH -n $JOBS_PER_CHUNK
#SBATCH -A hpc_csbg_25

# Load necessary modules/environments
module load intel/2021.5.0
module load intel-mpi/2021.5.1
module load parallel/20220522/intel-2021.5.0
source activate biopython  # Make sure this environment exists and has Biopython

mkdir -p $OUTPUT_DIR # Ensure output directory exists

# The core parallel command:
# -a \$chunk reads the list of input files from the chunk file.
# The command uses shell substitution to pass the filename into the embedded Python script.
parallel --wd $SLRUM_RUN_DIR \
  --jobs $JOBS_PER_CHUNK -a $chunk \
  '
    infile="{}"
    base_name=\$(basename "\$infile" .faa)
    outfile="${OUTPUT_DIR}/\$base_name.faa"
    
    python - <<PY_EOF
from Bio import SeqIO
infile = "\$infile"
outfile = "\$outfile"
with open(outfile, "w") as out:
    SeqIO.write(
        (r for r in SeqIO.parse(infile, "fasta") if len(r.seq) > 50),
        out,
        "fasta"
    )
PY_EOF
'

date
exit 0

EOF

  chmod +x "$SCRIPT_DIR/${base}.sh"
  echo "Created Slurm script: $SCRIPT_DIR/${base}.sh"

done

echo "Done creating Slurm job scripts in $SCRIPT_DIR."
echo "Use 'sbatch $SCRIPT_DIR/part_00000.sh' to submit the first job."
