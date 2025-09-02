# Author: Gopal Srivastava
# Date: Sept 2, 2025

#!/bin/bash

# === Configuration ===
IDS="lst"      # master input file
CHUNK_DIR="./chunks"   # where chunk files will live
SCRIPT_DIR="./slurm_scripts"  # where per-chunk Slurm scripts will live
LINES_PER_CHUNK=2000      # adjust to control number of chunks
JOBS_PER_CHUNK=20         # number of concurrent tasks per chunk

SLRUM_RUN_DIR=$( pwd )

mkdir -p "$CHUNK_DIR" "$SCRIPT_DIR"

# Clean old chunk + scripts
rm -f "$CHUNK_DIR"/part_* "$SCRIPT_DIR"/chunk_*.sh

# === Split into chunks ===
split -l "$LINES_PER_CHUNK" --numeric-suffixes=0 --suffix-length=5 \
  "$IDS" "$CHUNK_DIR/part_"

# === Make a Slurm script for each chunk ===
for chunk in "$CHUNK_DIR"/part_*; do
  base=$(basename "$chunk")
  jobname="bvbrc_${base}"

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

module load intel/2021.5.0
module load intel-mpi/2021.5.1
module load parallel/20220522/intel-2021.5.0


mkdir -p faa_cleaned

parallel --wd $SLRUM_RUN_DIR \
 --jobs $JOBS_PER_CHUNK -a $chunk \
 'python RemoveUndefinedorHypotheticalProteins.py -f {}'

date
exit 0

EOF

chmod +x "$SCRIPT_DIR/${base}.sh"

done
