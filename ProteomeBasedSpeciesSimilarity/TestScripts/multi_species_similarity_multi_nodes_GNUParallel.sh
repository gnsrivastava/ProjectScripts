#!/bin/bash 
: <<'END_COMMENT'
***********************************
Author: Gopal Srivastava
Date: Sept 8, 2025
************************************
Steps:

Prepare a list of species (fasta files) in a directory.
For each species, create a Diamond database.
Then for each pair (i, j) where i != j, run the two Diamond searches and compute the average.

# ================================================
Pre-requisites:
Install Diamond for faster blastp (https://github.com/bbuchfink/diamond)
Install Python3

# ============================================
Steps to get species similarity using diamond 

1. Precompute databases for each species.
For each unordered pair (i, j) with i != j, do:
a. Run Diamond: i vs j -> get top hits for each i, average the identities -> avg_i_j
b. Run Diamond: j vs i -> get top hits for each j, average the identities -> avg_j_i
c. Overall similarity for pair (i,j) = (avg_i_j + avg_j_i) / 2
Store the result in a matrix.

Check line number 50, 70 & 71 to change delimiter according to your input filename
END_COMMENT

# Check if correct number of arguments are provided
if [ $# -lt 2 ]; then
    echo "Usage: $0 <input_directory> <output_prefix> [num_chunks]"
    echo "Input directory should contain FASTA files for each species"
    exit 1
fi

INPUT_DIR=$1
PREFIX=$2
NUM_CHUNKS=${3:-100}  # Default to 100 chunks if not specified

# Create output directory
OUTPUT_DIR="${PREFIX}_output"
mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/databases"
mkdir -p "$OUTPUT_DIR/results"
mkdir -p "$OUTPUT_DIR/chunks"
mkdir -p "$OUTPUT_DIR/slurm_scripts"
mkdir -p "$OUTPUT_DIR/slurm_logs"

# Get list of all FASTA files
FASTA_FILES=($(find "$INPUT_DIR" -name "*.faa" -o -name "*.fasta" -o -name "*.fa"))
if [ ${#FASTA_FILES[@]} -lt 2 ]; then
    echo "Error: Need at least 2 FASTA files for comparison"
    exit 1
fi

echo "Found ${#FASTA_FILES[@]} FASTA files"

# Create Diamond databases for all files
echo "Creating Diamond databases..."
for file in "${FASTA_FILES[@]}"; do
    base_name=$(basename "$file" | cut -d'_' -f1) # My files are named as genomeid_nonredundant.faa. Replace delimiter according to the filename
    db_file="$OUTPUT_DIR/databases/${base_name}.dmnd"
    
    if [ ! -f "$db_file" ]; then
        echo "Creating database for $base_name..."
        diamond makedb --in "$file" -d "$db_file" --quiet
    fi
done

# Generate all possible pairs of species
echo "Generating species pairs..."
PAIRS_FILE="$OUTPUT_DIR/species_pairs.txt"
> "$PAIRS_FILE"  # Clear the file

for i in "${!FASTA_FILES[@]}"; do
    for j in "${!FASTA_FILES[@]}"; do
        if [ $i -lt $j ]; then  # Only process each pair once
            file1="${FASTA_FILES[$i]}"
            file2="${FASTA_FILES[$j]}"
            base1=$(basename "$file1" | cut -d'_' -f1) # My files are named as genomeid_nonredundant.faa. Replace delimiter according to the filename
            base2=$(basename "$file2" | cut -d'_' -f1) # My files are named as genomeid_nonredundant.faa. Replace delimiter according to the filename
            echo "$file1 $file2 $base1 $base2" >> "$PAIRS_FILE"
        fi
    done
done

# Split the pairs file into chunks
echo "Splitting pairs file into $NUM_CHUNKS chunks..."
TOTAL_PAIRS=$(wc -l < "$PAIRS_FILE")
PAIRS_PER_CHUNK=$(( (TOTAL_PAIRS + NUM_CHUNKS - 1) / NUM_CHUNKS ))  # Ceiling division

split -l $PAIRS_PER_CHUNK --numeric-suffixes=1 --suffix-length=3 \
    "$PAIRS_FILE" "$OUTPUT_DIR/chunks/chunk_"

# Create a SLURM script template
SLURM_SCRIPT_TEMPLATE="$OUTPUT_DIR/slurm_script_template.sh"
cat > "$SLURM_SCRIPT_TEMPLATE" << 'EOF'
#!/bin/bash
#SBATCH -J species_similarity_chunk_${CHUNK_NUM}
#SBATCH -o ${OUTPUT_DIR}/slurm_logs/slurm_%j.out
#SBATCH -e ${OUTPUT_DIR}/slurm_logs/slurm_%j.err
#SBATCH -p workq
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -A hpc_csbg23

module load intel/2021.5.0
module load intel-mpi/2021.5.1
module load diamond/2.0.14

# Function to process a single pair
process_pair() {
    file1=$1
    file2=$2
    base1=$3
    base2=$4
    
    # Create temporary directory for this pair
    TMPDIR=$(mktemp -d)
    
    # Run Diamond in both directions
    diamond blastp -q "$file1" -d "${OUTPUT_DIR}/databases/${base2}.dmnd" -o "$TMPDIR/${base1}_vs_${base2}.tsv" \
        --id -f 6 qseqid sseqid pident -k 1 --quiet --fast
    
    diamond blastp -q "$file2" -d "${OUTPUT_DIR}/databases/${base1}.dmnd" -o "$TMPDIR/${base2}_vs_${base1}.tsv" \
        --id -f 6 qseqid sseqid pident -k 1 --quiet --fast
    
    # Calculate average similarities
    avg_AB=$(awk '{sum += $3; count++} END {if (count>0) print sum/count; else print 0}' "$TMPDIR/${base1}_vs_${base2}.tsv")
    avg_BA=$(awk '{sum += $3; count++} END {if (count>0) print sum/count; else print 0}' "$TMPDIR/${base2}_vs_${base1}.tsv")
    
    # Calculate final average
    final_avg=$(echo "scale=4; ($avg_AB + $avg_BA) / 2" | bc)
    
    # Count sequences in each species
    countA=$(grep -c ">" "$file1")
    countB=$(grep -c ">" "$file2")
    
    # Save results
    echo "$base1,$base2,$countA,$countB,$avg_AB,$avg_BA,$final_avg" >> "${OUTPUT_DIR}/results/pairwise_results_chunk_${CHUNK_NUM}.csv"
    
    # Clean up
    rm -rf "$TMPDIR"
}

export -f process_pair
export OUTPUT_DIR

# Process all pairs in the chunk
echo "Processing chunk ${CHUNK_NUM} with $(wc -l < ${CHUNK_FILE}) pairs"
cat "${CHUNK_FILE}" | parallel -j 20 --colsep ' ' 'process_pair {1} {2} {3} {4}'

echo "Chunk ${CHUNK_NUM} processing complete"
EOF

# Create SLURM scripts for each chunk
echo "Creating SLURM scripts for each chunk..."
for chunk_file in "$OUTPUT_DIR/chunks/chunk_"*; do
    chunk_num=$(echo "$chunk_file" | grep -o '[0-9][0-9][0-9]$')
    
    # Create the SLURM script for this chunk
    slurm_script="$OUTPUT_DIR/slurm_scripts/process_chunk_${chunk_num}.sh"
    
    # Replace placeholders in the template
    sed "s/\${CHUNK_NUM}/${chunk_num}/g; s|\${OUTPUT_DIR}|${OUTPUT_DIR}|g; s|\${CHUNK_FILE}|${chunk_file}|g" \
        "$SLURM_SCRIPT_TEMPLATE" > "$slurm_script"
    
    chmod +x "$slurm_script"
    echo "Created SLURM script: $slurm_script"
done

# Create a master submission script
MASTER_SUBMIT_SCRIPT="$OUTPUT_DIR/submit_all_chunks.sh"
cat > "$MASTER_SUBMIT_SCRIPT" << EOF
#!/bin/bash
# Master script to submit all chunk processing jobs

for script in "$OUTPUT_DIR/slurm_scripts/process_chunk_"*.sh; do
    echo "Submitting \$script"
    sbatch "\$script"
    sleep 1
done

echo "All chunk jobs submitted. Use 'squeue -u \$USER' to check status."
EOF

chmod +x "$MASTER_SUBMIT_SCRIPT"

# Create a script to merge results after all jobs complete
MERGE_SCRIPT="$OUTPUT_DIR/merge_results.sh"
cat > "$MERGE_SCRIPT" << 'EOF'
#!/bin/bash
# Script to merge results from all chunks

# Check if all chunk results are available
expected_chunks=$(ls -1 "${OUTPUT_DIR}/chunks/chunk_"* | wc -l)
available_results=$(ls -1 "${OUTPUT_DIR}/results/pairwise_results_chunk_"*.csv 2>/dev/null | wc -l)

if [ "$available_results" -lt "$expected_chunks" ]; then
    echo "Warning: Only $available_results out of $expected_chunks chunks have completed."
    echo "Some jobs may still be running or failed."
    read -p "Continue with merging anyway? (y/N) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        exit 1
    fi
fi

# Create header for the combined results file
echo "SpeciesA,SpeciesB,CountA,CountB,Avg_A_to_B,Avg_B_to_A,Final_Avg" > "${OUTPUT_DIR}/results/pairwise_results_combined.csv"

# Combine all chunk results (skip headers)
for result_file in "${OUTPUT_DIR}/results/pairwise_results_chunk_"*.csv; do
    # Skip the header line in each file
    tail -n +2 "$result_file" >> "${OUTPUT_DIR}/results/pairwise_results_combined.csv"
done

# Generate similarity matrix
echo "Generating similarity matrix..."
python3 -c "
import pandas as pd
import numpy as np

# Read results
df = pd.read_csv('${OUTPUT_DIR}/results/pairwise_results_combined.csv')

# Get unique species
species = sorted(set(df['SpeciesA'].unique()) | set(df['SpeciesB'].unique()))

# Create empty matrix
matrix = pd.DataFrame(np.zeros((len(species), len(species))), index=species, columns=species)

# Fill matrix with similarity values
for _, row in df.iterrows():
    spA, spB, avg = row['SpeciesA'], row['SpeciesB'], row['Final_Avg']
    matrix.loc[spA, spB] = avg
    matrix.loc[spB, spA] = avg

# Set diagonal to 100 (self-similarity)
np.fill_diagonal(matrix.values, 100)

# Save matrix
matrix.to_csv('${OUTPUT_DIR}/results/similarity_matrix.csv')
print('Similarity matrix saved to ${OUTPUT_DIR}/results/similarity_matrix.csv')
"

echo "Results merged successfully!"
echo "Combined results: ${OUTPUT_DIR}/results/pairwise_results_combined.csv"
echo "Similarity matrix: ${OUTPUT_DIR}/results/similarity_matrix.csv"
EOF

chmod +x "$MERGE_SCRIPT"

echo "Setup complete!"
echo "1. Submit all chunk jobs: ${MASTER_SUBMIT_SCRIPT}"
echo "2. After all jobs complete, merge results: ${MERGE_SCRIPT}"
echo "3. Check SLURM logs in: ${OUTPUT_DIR}/slurm_logs/"
