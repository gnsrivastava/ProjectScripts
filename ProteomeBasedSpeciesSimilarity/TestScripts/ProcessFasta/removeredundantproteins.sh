#!/bin/bash

# === Configuration ===
IDS="lst"                      # master input file (list of FASTA files)
CHUNK_DIR="./chunks"           # where chunk files will live
SCRIPT_DIR="./slurm_scripts"   # where per-chunk Slurm scripts will live
OUTPUT_DIR="./removed_redundant" # output directory for processed files
LINES_PER_CHUNK=2000           # adjust to control number of chunks
JOBS_PER_CHUNK=20              # number of concurrent tasks per chunk
NUM_JOBS=20                    # number of parallel CD-HIT jobs

SLRUM_RUN_DIR=$(pwd)

# Create necessary directories
mkdir -p "$CHUNK_DIR" "$SCRIPT_DIR" "$OUTPUT_DIR"

# Clean old chunk + scripts
rm -f "$CHUNK_DIR"/part_* "$SCRIPT_DIR"/chunk_*.sh

# === Split into chunks ===
echo "Splitting master file into chunks..."
split -l "$LINES_PER_CHUNK" --numeric-suffixes=0 --suffix-length=5 \
  "$IDS" "$CHUNK_DIR/part_"

# === Make a Slurm script for each chunk ===
echo "Creating Slurm scripts for each chunk..."
for chunk in "$CHUNK_DIR"/part_*; do
  base=$(basename "$chunk")
  jobname="cdhit_${base}"

  cat > "$SCRIPT_DIR/${base}.sh" <<EOF
#!/bin/bash

#SBATCH -J $jobname
#SBATCH -o $jobname.%j.out
#SBATCH -e $jobname.%j.err
#SBATCH -p workq
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n $JOBS_PER_CHUNK
#SBATCH -A hpc_csbg23

source activate cdhit

module load intel/2021.5.0
module load intel-mpi/2021.5.1
module load parallel/20220522/intel-2021.5.0

# Process all FASTA files in the chunk in parallel
parallel --wd $SLRUM_RUN_DIR -j $NUM_JOBS -a $chunk '
file={};
base_name=\$(basename "\$file" .faa);
input_file="./filtered/\$file";
output_file="$OUTPUT_DIR/\${base_name}_nonredundant.faa";

# Count sequences in the input file
seq_count=\$(grep -c ">" "\$input_file");

if [ "\$seq_count" -gt 1 ]; then
    echo "Processing \$input_file with \$seq_count sequences...";
    cd-hit -i "\$input_file" -o "\$output_file" -c 0.9 -n 5 -M 16000 -d 0;
    echo "Completed processing \$input_file";
else
    echo "Copying \$input_file (only \$seq_count sequences)...";
    cp "\$input_file" "\$output_file";
    echo "Copied \$input_file to \$output_file";
fi
'

date
exit 0
EOF

  chmod +x "$SCRIPT_DIR/${base}.sh"
  echo "Created Slurm script: $SCRIPT_DIR/${base}.sh"
done

echo "All Slurm scripts have been created. Submit them using:"
echo "for script in $SCRIPT_DIR/part_*.sh; do sbatch \$script; sleep 10; done"
