# Usage of the script: parallel -j 4 python calculate_dipeptide_frequencies.py {} ../DIPEP_RES/{/.}.csv ::: *.fa
# Adjust the number of workers (j) according to availability. 
import pandas as pd
from Bio import SeqIO
from itertools import product
from collections import Counter
import argparse
import os

def calculate_dipeptide_frequencies(fasta_file, output_csv):
    """
    Calculates frequencies for all 400 possible dipeptides from a FASTA file.
    Formula: Frequency = (Count of Dipeptide) / (Length of Protein)
    """
    # 20 standard amino acids
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    
    # Generate all 400 possible dipeptides (AA, AC, AD...)
    all_dipeptides = [''.join(p) for p in product(amino_acids, repeat=2)]
    
    results = []

    if not os.path.exists(fasta_file):
        print(f"Error: FASTA file not found: {fasta_file}")
        return

    print(f"Reading sequences from {fasta_file}...")
    
    # Parse FASTA and calculate frequencies
    count_processed = 0
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq).upper()
        seq_length = len(sequence)
        
        # A dipeptide requires at least 2 amino acids
        if seq_length < 2:
            print(f"Skipping {record.id}: Sequence length ({seq_length}) too short.")
            continue

        # Extract all overlapping dipeptides
        # Range is length-1 because the last amino acid cannot start a dipeptide
        dipeps_in_seq = [sequence[i:i+2] for i in range(seq_length - 1)]
        counts = Counter(dipeps_in_seq)
        
        # Calculate frequency for each of the 400 possible dipeptides
        # Formula: count / length
        row = {'Protein_ID': record.id, 'Sequence_Length': seq_length}
        for dp in all_dipeptides:
            # Using .get(dp, 0) returns 0 if the dipeptide was not found in the sequence
            row[dp] = counts.get(dp, 0) / seq_length
            
        results.append(row)
        count_processed += 1

    # Save to CSV
    if results:
        df = pd.DataFrame(results)
        df.to_csv(output_csv, index=False)
        print(f"Successfully processed {count_processed} proteins.")
        print(f"Results saved to: {output_csv}")
    else:
        print("No valid sequences were found in the input file.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate 400 dipeptide frequencies from a FASTA file.")
    parser.add_argument("input", help="Path to the input FASTA file")
    parser.add_argument("output", help="Path to the output CSV file")
    
    args = parser.parse_args()
    
    calculate_dipeptide_frequencies(args.input, args.output)
