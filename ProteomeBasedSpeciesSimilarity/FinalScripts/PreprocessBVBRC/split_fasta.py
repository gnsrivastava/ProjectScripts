"""
make sure to activate biopython conda environment or have biopython installed
Split a FASTA file into two parts.

Multiple splitting options:
1. By number of sequences (first N sequences vs rest)
2. By ratio (e.g., 50/50, 80/20)
3. By random split (for train/test sets)
4. By sequence IDs (from a list)

Requirements:
    pip install biopython
"""

import os
import random
from pathlib import Path
from typing import List, Tuple, Optional
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse


def split_by_count(
    input_file: str,
    n_first: int,
    output_prefix: str = None,
    verbose: bool = True
) -> Tuple[str, str]:
    """
    Split FASTA file: first N sequences vs remaining sequences.
    
    Args:
        input_file: Input FASTA file
        n_first: Number of sequences in first part
        output_prefix: Prefix for output files (default: input filename)
        verbose: Print progress
    
    Returns:
        Tuple of (part1_file, part2_file)
    """
    if output_prefix is None:
        output_prefix = Path(input_file).stem
    
    suffix = Path(input_file).suffix or '.fa'
    part1_file = f"{output_prefix}_part1{suffix}"
    part2_file = f"{output_prefix}_part2{suffix}"
    
    records = list(SeqIO.parse(input_file, 'fasta'))
    total = len(records)
    
    if n_first >= total:
        raise ValueError(f"n_first ({n_first}) >= total sequences ({total})")
    
    part1 = records[:n_first]
    part2 = records[n_first:]
    
    SeqIO.write(part1, part1_file, 'fasta')
    SeqIO.write(part2, part2_file, 'fasta')
    
    if verbose:
        print(f"Split by count: {n_first} / {total - n_first}")
        print(f"  Part 1: {part1_file} ({len(part1)} sequences)")
        print(f"  Part 2: {part2_file} ({len(part2)} sequences)")
    
    return part1_file, part2_file


def split_by_ratio(
    input_file: str,
    ratio: float = 0.5,
    output_prefix: str = None,
    verbose: bool = True
) -> Tuple[str, str]:
    """
    Split FASTA file by ratio (e.g., 0.5 = 50/50 split).
    
    Args:
        input_file: Input FASTA file
        ratio: Fraction of sequences in first part (0.0 to 1.0)
        output_prefix: Prefix for output files
        verbose: Print progress
    
    Returns:
        Tuple of (part1_file, part2_file)
    """
    if not 0 < ratio < 1:
        raise ValueError(f"Ratio must be between 0 and 1, got {ratio}")
    
    records = list(SeqIO.parse(input_file, 'fasta'))
    total = len(records)
    n_first = int(total * ratio)
    
    if n_first == 0:
        n_first = 1
    if n_first == total:
        n_first = total - 1
    
    return split_by_count(input_file, n_first, output_prefix, verbose)


def split_random(
    input_file: str,
    ratio: float = 0.5,
    output_prefix: str = None,
    seed: int = None,
    verbose: bool = True
) -> Tuple[str, str]:
    """
    Randomly split FASTA file (useful for train/test splits).
    
    Args:
        input_file: Input FASTA file
        ratio: Fraction of sequences in first part
        output_prefix: Prefix for output files
        seed: Random seed for reproducibility
        verbose: Print progress
    
    Returns:
        Tuple of (part1_file, part2_file)
    """
    if seed is not None:
        random.seed(seed)
    
    if output_prefix is None:
        output_prefix = Path(input_file).stem
    
    suffix = Path(input_file).suffix or '.fa'
    part1_file = f"{output_prefix}_part1{suffix}"
    part2_file = f"{output_prefix}_part2{suffix}"
    
    records = list(SeqIO.parse(input_file, 'fasta'))
    total = len(records)
    
    # Shuffle indices
    indices = list(range(total))
    random.shuffle(indices)
    
    n_first = int(total * ratio)
    if n_first == 0:
        n_first = 1
    if n_first == total:
        n_first = total - 1
    
    part1_indices = set(indices[:n_first])
    
    part1 = [records[i] for i in range(total) if i in part1_indices]
    part2 = [records[i] for i in range(total) if i not in part1_indices]
    
    SeqIO.write(part1, part1_file, 'fasta')
    SeqIO.write(part2, part2_file, 'fasta')
    
    if verbose:
        print(f"Random split (ratio={ratio}, seed={seed}): {len(part1)} / {len(part2)}")
        print(f"  Part 1: {part1_file} ({len(part1)} sequences)")
        print(f"  Part 2: {part2_file} ({len(part2)} sequences)")
    
    return part1_file, part2_file


def split_by_ids(
    input_file: str,
    ids_part1: List[str],
    output_prefix: str = None,
    verbose: bool = True
) -> Tuple[str, str]:
    """
    Split FASTA file by sequence IDs.
    
    Args:
        input_file: Input FASTA file
        ids_part1: List of sequence IDs for part 1
        output_prefix: Prefix for output files
        verbose: Print progress
    
    Returns:
        Tuple of (part1_file, part2_file)
    """
    if output_prefix is None:
        output_prefix = Path(input_file).stem
    
    suffix = Path(input_file).suffix or '.fa'
    part1_file = f"{output_prefix}_part1{suffix}"
    part2_file = f"{output_prefix}_part2{suffix}"
    
    ids_set = set(ids_part1)
    
    part1 = []
    part2 = []
    
    for record in SeqIO.parse(input_file, 'fasta'):
        if record.id in ids_set:
            part1.append(record)
        else:
            part2.append(record)
    
    SeqIO.write(part1, part1_file, 'fasta')
    SeqIO.write(part2, part2_file, 'fasta')
    
    if verbose:
        print(f"Split by IDs: {len(part1)} / {len(part2)}")
        print(f"  Part 1: {part1_file} ({len(part1)} sequences)")
        print(f"  Part 2: {part2_file} ({len(part2)} sequences)")
        if len(part1) != len(ids_part1):
            print(f"  Warning: {len(ids_part1) - len(part1)} IDs not found in input file")
    
    return part1_file, part2_file


def split_by_ids_file(
    input_file: str,
    ids_file: str,
    output_prefix: str = None,
    verbose: bool = True
) -> Tuple[str, str]:
    """
    Split FASTA file using IDs from a file (one ID per line).
    
    Args:
        input_file: Input FASTA file
        ids_file: File with sequence IDs (one per line)
        output_prefix: Prefix for output files
        verbose: Print progress
    
    Returns:
        Tuple of (part1_file, part2_file)
    """
    with open(ids_file) as f:
        ids = [line.strip() for line in f if line.strip()]
    
    return split_by_ids(input_file, ids, output_prefix, verbose)


def split_half(
    input_file: str,
    output_prefix: str = None,
    verbose: bool = True
) -> Tuple[str, str]:
    """
    Split FASTA file into two equal halves.
    
    Args:
        input_file: Input FASTA file
        output_prefix: Prefix for output files
        verbose: Print progress
    
    Returns:
        Tuple of (part1_file, part2_file)
    """
    return split_by_ratio(input_file, 0.5, output_prefix, verbose)


def split_into_n_parts(
    input_file: str,
    n_parts: int,
    output_prefix: str = None,
    verbose: bool = True
) -> List[str]:
    """
    Split FASTA file into N equal parts.
    
    Args:
        input_file: Input FASTA file
        n_parts: Number of parts to split into
        output_prefix: Prefix for output files
        verbose: Print progress
    
    Returns:
        List of output file paths
    """
    if output_prefix is None:
        output_prefix = Path(input_file).stem
    
    suffix = Path(input_file).suffix or '.fa'
    
    records = list(SeqIO.parse(input_file, 'fasta'))
    total = len(records)
    
    if n_parts > total:
        raise ValueError(f"n_parts ({n_parts}) > total sequences ({total})")
    
    # Calculate sequences per part
    base_size = total // n_parts
    remainder = total % n_parts
    
    output_files = []
    start = 0
    
    for i in range(n_parts):
        # Distribute remainder across first few parts
        part_size = base_size + (1 if i < remainder else 0)
        end = start + part_size
        
        part_records = records[start:end]
        output_file = f"{output_prefix}_part{i+1}{suffix}"
        
        SeqIO.write(part_records, output_file, 'fasta')
        output_files.append(output_file)
        
        if verbose:
            print(f"  Part {i+1}: {output_file} ({len(part_records)} sequences)")
        
        start = end
    
    if verbose:
        print(f"\nSplit into {n_parts} parts from {total} sequences")
    
    return output_files


# =============================================================================
# COMMAND LINE INTERFACE
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Split a FASTA file into two (or more) parts',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Split in half (50/50)
  python split_fasta.py input.fa
  
  # Split first 100 sequences vs rest
  python split_fasta.py input.fa -n 100
  
  # Split by ratio (80/20)
  python split_fasta.py input.fa -r 0.8
  
  # Random split for train/test (70/30)
  python split_fasta.py input.fa -r 0.7 --random --seed 42
  
  # Split by IDs from file
  python split_fasta.py input.fa --ids-file selected_ids.txt
  
  # Split into 5 equal parts
  python split_fasta.py input.fa --n-parts 5
  
  # Custom output prefix
  python split_fasta.py input.fa -o my_split
        """
    )
    
    parser.add_argument('input', help='Input FASTA file')
    parser.add_argument('-o', '--output', help='Output prefix (default: input filename)')
    parser.add_argument('-n', '--count', type=int, 
                        help='Number of sequences in first part')
    parser.add_argument('-r', '--ratio', type=float, default=0.5,
                        help='Ratio for first part (default: 0.5)')
    parser.add_argument('--random', action='store_true',
                        help='Random split (for train/test)')
    parser.add_argument('--seed', type=int, 
                        help='Random seed for reproducibility')
    parser.add_argument('--ids-file', 
                        help='File with sequence IDs for part 1 (one per line)')
    parser.add_argument('--n-parts', type=int,
                        help='Split into N equal parts')
    parser.add_argument('-q', '--quiet', action='store_true',
                        help='Suppress output messages')
    
    args = parser.parse_args()
    verbose = not args.quiet
    
    if args.n_parts:
        # Split into N parts
        split_into_n_parts(args.input, args.n_parts, args.output, verbose)
    
    elif args.ids_file:
        # Split by IDs from file
        split_by_ids_file(args.input, args.ids_file, args.output, verbose)
    
    elif args.count:
        # Split by count
        split_by_count(args.input, args.count, args.output, verbose)
    
    elif args.random:
        # Random split
        split_random(args.input, args.ratio, args.output, args.seed, verbose)
    
    else:
        # Split by ratio (default: 50/50)
        split_by_ratio(args.input, args.ratio, args.output, verbose)


# =============================================================================
# EXAMPLE USAGE
# =============================================================================

if __name__ == "__main__":
    import sys
    
    if len(sys.argv) > 1:
        main()
    else:
        # Demo
        print("=" * 60)
        print("FASTA Splitter")
        print("=" * 60)
        print("\nUsage examples:")
        print("  python split_fasta.py input.fa                    # 50/50 split")
        print("  python split_fasta.py input.fa -n 100             # First 100 vs rest")
        print("  python split_fasta.py input.fa -r 0.8             # 80/20 split")
        print("  python split_fasta.py input.fa -r 0.7 --random    # Random 70/30")
        print("  python split_fasta.py input.fa --n-parts 5        # Split into 5 parts")
        print("  python split_fasta.py input.fa --ids-file ids.txt # Split by ID list")
        print("\nRun with -h for full help.")
