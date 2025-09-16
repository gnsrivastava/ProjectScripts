#!/bin/python3
# This script returns species vs species similarity matrix with diagona representing similarity between same species

import pandas as pd
from pathlib import Path

# Data path
DATA_DIR = Path(".")

# Load number of sequences per species
numberofseqs = pd.read_csv('../numberofseqs.txt', sep=':')
numberofseqs.Species = numberofseqs.Species.astype(str)

# Dictionary to look up number of sequences by species ID
seq_count = dict(zip(numberofseqs.Species, numberofseqs.NumberofSeqs))

# Store results
results = {}

# Iterate over each hungarian output file (*.csv)
for filename in DATA_DIR.glob("*.csv"):
    fname = filename.name

    # Extract species IDs
    if "_vs_" not in fname:
        continue
    a, b = fname.split("_vs_")
    b = b.replace(".tsv_hungarian.csv", "").replace(".tsv", "")

    # If same species → force similarity = 100
    if a == b:
        results.setdefault(a, {})[b] = 100.0
        continue

    # Get sequence counts
    if a not in seq_count or b not in seq_count:
        print(f"Skipping {fname}: missing seq count for {a} or {b}")
        continue
    na, nb = seq_count[a], seq_count[b]

    # Minimum number of sequences
    min_matches = min(na, nb)

    # Read the file and compute similarity
    df = pd.read_csv(filename)
    if "avg_pident" not in df.columns:
        print(f"Skipping {fname}: no pident column")
        continue

    # If df has fewer alignments than expected, pad with zeros
    observed = len(df)
    if observed < min_matches:
        missing = min_matches - observed
        speciesSimilarity = (df["avg_pident"].sum() + missing * df.avg_pident.min()) / min_matches
    else:
        speciesSimilarity = df["avg_pident"].sum() / min_matches

    # Store in nested dict
    results.setdefault(a, {})[b] = speciesSimilarity

# Convert results to DataFrame (rows = A, cols = B)
similarity_matrix = pd.DataFrame.from_dict(results, orient="index")

# --- Symmetrize: fill NaNs with transpose values ---
similarity_matrix = similarity_matrix.combine_first(similarity_matrix.T)

# Ensure diagonal = 100 (in case no self files existed)
for sp in similarity_matrix.index:
    similarity_matrix.loc[sp, sp] = 100.0

# Save matrix
similarity_matrix.to_csv("species_similarity_matrix.csv")

print("Wrote species_similarity_matrix.csv")
