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


# Check if correct number of arguments are provided
if [ $# -lt 2 ]; then
    echo "Usage: $0 <input_directory> <output_prefix> [num_parallel_jobs]"
    echo "Input directory should contain FASTA files for each species"
    exit 1
fi

INPUT_DIR=$1
PREFIX=$2
NUM_JOBS=${3:-20}  # Default to 20 parallel jobs if not specified

# Create output directory
OUTPUT_DIR="${PREFIX}_output"
mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/databases"
mkdir -p "$OUTPUT_DIR/results"

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
    base_name=$(basename "$file" | cut -d. -f1)
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
            base1=$(basename "$file1" | cut -d. -f1)
            base2=$(basename "$file2" | cut -d. -f1)
            echo "$file1 $file2 $base1 $base2" >> "$PAIRS_FILE"
        fi
    done
done

# Function to process a single pair
process_pair() {
    file1=$1
    file2=$2
    base1=$3
    base2=$4
    
    # Create temporary directory for this pair
    TMPDIR=$(mktemp -d)
    
    # Run Diamond in both directions
    diamond blastp -q "$file1" -d "$OUTPUT_DIR/databases/${base2}.dmnd" -o "$TMPDIR/${base1}_vs_${base2}.tsv" \
        --id -f 6 qseqid sseqid pident -k 1 --quiet --fast
    
    diamond blastp -q "$file2" -d "$OUTPUT_DIR/databases/${base1}.dmnd" -o "$TMPDIR/${base2}_vs_${base1}.tsv" \
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
    echo "$base1,$base2,$countA,$countB,$avg_AB,$avg_BA,$final_avg" >> "$OUTPUT_DIR/results/pairwise_results.csv"
    
    # Copy result files if needed
    cp "$TMPDIR/${base1}_vs_${base2}.tsv" "$OUTPUT_DIR/results/"
    cp "$TMPDIR/${base2}_vs_${base1}.tsv" "$OUTPUT_DIR/results/"
    
    # Clean up
    rm -rf "$TMPDIR"
    
    echo "Processed pair: $base1 vs $base2"
}

export -f process_pair
export OUTPUT_DIR

# Create results file header
echo "SpeciesA,SpeciesB,CountA,CountB,Avg_A_to_B,Avg_B_to_A,Final_Avg" > "$OUTPUT_DIR/results/pairwise_results.csv"

# Process all pairs in parallel
echo "Processing all species pairs..."
cat "$PAIRS_FILE" | parallel -j "$NUM_JOBS" --colsep ' ' 'process_pair {1} {2} {3} {4}'

# Generate similarity matrix
echo "Generating similarity matrix..."
python -c "
import pandas as pd
import numpy as np

# Read results
df = pd.read_csv('$OUTPUT_DIR/results/pairwise_results.csv')

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
matrix.to_csv('$OUTPUT_DIR/results/similarity_matrix.csv')
print('Similarity matrix saved to $OUTPUT_DIR/results/similarity_matrix.csv')
"

echo "Analysis complete. Results saved to $OUTPUT_DIR/results/"
echo "Pairwise results: $OUTPUT_DIR/results/pairwise_results.csv"
echo "Similarity matrix: $OUTPUT_DIR/results/similarity_matrix.csv"
