#!/usr/bin/env python3
# mpi_tanimoto_allvsall.py
#
# All-vs-all Tanimoto with OpenMPI + mpi4py + RDKit (batched rows x batched cols).
#
# Usage (32 ranks):
#   mpirun -np 32 python mpi_tanimoto_allvsall.py \
#       --in smiles.tsv \
#       --out-matrix tanimoto_matrix.csv \
#       --out-pairs tanimoto_pairs.tsv \
#       --pairs-thresh 0.85 \
#       --row-batch 2000 --col-batch 2000
#
# Input format: TSV/CSV with columns: name, canonical_smiles (or SMILES)
#
# Notes:
# - Fingerprint: Morgan radius=2, nBits=2048
# - Invalid SMILES -> NaN in matrix
# - Work splitting: round-robin over row batches across ranks
# - The pairs TSV is written only on rank 0 after gathering

from __future__ import annotations
import argparse
import math
import os
from typing import List, Tuple

import numpy as np
import pandas as pd
from mpi4py import MPI

from rdkit import Chem, DataStructs
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator


# ---------- MPI setup ----------
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


# ---------- Args ----------
def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--in", dest="infile", required=True,
                   help="Input TSV/CSV with columns: name, canonical_smiles (or SMILES)")
    p.add_argument("--out-matrix", required=True,
                   help="Output CSV of symmetric Tanimoto matrix (names as header/index)")
    p.add_argument("--out-pairs", default=None,
                   help="Optional output TSV of (i,j,name_i,name_j,tanimoto)")
    p.add_argument("--pairs-thresh", type=float, default=0.85,
                   help="Minimum similarity to write to --out-pairs (default 0.85)")
    p.add_argument("--row-batch", type=int, default=2000,
                   help="Row batch size (default 2000)")
    p.add_argument("--col-batch", type=int, default=2000,
                   help="Column batch size (default 2000)")
    return p.parse_args()


# ---------- IO ----------
def read_table(path: str) -> pd.DataFrame:
    # Autodetect delimiter
    try:
        df = pd.read_csv(path)
    except Exception:
        df = pd.read_csv(path, sep="\t")
    cols = {c.strip().lower(): c for c in df.columns}
    #name_col = cols.get("name")
    #smi_col = cols.get("canonical_smiles") or cols.get("smiles")
    if not name_col or not smi_col:
        raise ValueError("Input must have columns 'name' and 'canonical_smiles' (or 'SMILES').")
    out = df[['name', 'canonical_smiles']].rename(columns={'name': "name", 'canonical_smiles': "smiles"}).copy()
    out["name"] = out["name"].astype(str)
    out["smiles"] = out["smiles"].astype(str)
    return out


# ---------- FP ----------
morgan_gen = GetMorganGenerator(radius=2, fpSize=2048)

def smiles_to_fp(s: str):
    if not isinstance(s, str) or not s.strip():
        return None
    try:
        m = Chem.MolFromSmiles(s)
        if m is None:
            return None
        return morgan_gen.GetFingerprint(m)
    except Exception:
        return None


# ---------- Batching ----------
def make_batches(n: int, batch_size: int) -> List[Tuple[int, int]]:
    """Return list of (start, end) half-open index ranges covering 0..n."""
    return [(i, min(i + batch_size, n)) for i in range(0, n, batch_size)]


def compute_block(qfps, qvalid, tfps, tvalid) -> np.ndarray:
    """
    Compute a (len(qfps) x len(tfps)) block. Places NaN where either side is invalid.
    Uses BulkTanimotoSimilarity on valid targets for speed.
    """
    rows = []
    # Pre-collect indices of valid targets
    t_idx_valid = [j for j, ok in enumerate(tvalid) if ok]
    tfps_valid  = [tfps[j] for j in t_idx_valid]

    for i, (qfp, ok_q) in enumerate(zip(qfps, qvalid)):
        if not ok_q:
            rows.append([np.nan] * len(tfps))
            continue
        sims_valid = DataStructs.BulkTanimotoSimilarity(qfp, tfps_valid) if tfps_valid else []
        # expand to full-length row with NaN where target invalid
        row = [np.nan] * len(tfps)
        for k, j in enumerate(t_idx_valid):
            row[j] = sims_valid[k]
        rows.append(row)
        if i % 50 == 0 and i > 0:
            print(f"Rank {rank}: processed {i} rows in current block.", flush=True)
    return np.asarray(rows, dtype=np.float32)


def main():
    args = parse_args()

    # Rank 0 loads and broadcasts
    if rank == 0:
        df = read_table(args.infile)
        names = df["name"].tolist()
        smiles = df["smiles"].tolist()
        n = len(names)
        print(f"[rank 0] Loaded {n} compounds.")
    else:
        names = None
        smiles = None
        n = None

    names = comm.bcast(names, root=0)
    smiles = comm.bcast(smiles, root=0)
    n = comm.bcast(len(names) if rank == 0 else None, root=0)

    # Precompute fingerprints once per rank (read-only reuse)
    fps = [smiles_to_fp(s) for s in smiles]
    valid = np.array([fp is not None for fp in fps], dtype=bool)

    # Batches for rows and columns (all ranks agree on these)
    row_batches = make_batches(n, args.row_batch)
    col_batches = make_batches(n, args.col_batch)

    # Assign row batches round-robin across ranks
    my_row_batches = row_batches[rank::size]
    print(f"[rank {rank}] Assigned {len(my_row_batches)} row batches.", flush=True)

    # For assembling on root, we will send (row_start, block_matrix) per processed block
    my_blocks = []
    for (rs, re) in my_row_batches:
        qfps   = fps[rs:re]
        qvalid = valid[rs:re]
        row_block_parts = []
        for (cs, ce) in col_batches:
            tfps   = fps[cs:ce]
            tvalid = valid[cs:ce]
            block = compute_block(qfps, qvalid, tfps, tvalid)
            row_block_parts.append(block)
            print(f"[rank {rank}] Block rows [{rs}:{re}) x cols [{cs}:{ce}) -> {block.shape}", flush=True)
        # concat horizontally to full-width block for these rows
        full_row_block = np.hstack(row_block_parts)
        my_blocks.append((rs, full_row_block))

    # Gather all blocks on root
    gathered = comm.gather(my_blocks, root=0)

    if rank == 0:
        # Preallocate final matrix
        M = np.full((n, n), np.nan, dtype=np.float32)
        # Place each gathered block into the right rows
        for rank_blocks in gathered:
            for (rs, block) in rank_blocks:
                re = rs + block.shape[0]
                M[rs:re, :] = block

        # Enforce symmetry by averaging upper/lower where both defined; prefer max otherwise
        iu = np.triu_indices(n, 1)
        upper = M[iu]
        lower = M[(iu[1], iu[0])]
        # When both are finite, average; when one is NaN, take the other
        both = ~np.isnan(upper) & ~np.isnan(lower)
        only_u = ~np.isnan(upper) & np.isnan(lower)
        only_l = np.isnan(upper) & ~np.isnan(lower)

        sym_u = np.copy(upper)
        sym_l = np.copy(lower)
        sym_u[both] = (upper[both] + lower[both]) / 2.0
        sym_l[both] = sym_u[both]
        sym_l[only_u] = upper[only_u]
        sym_u[only_l] = lower[only_l]

        M[iu] = sym_u
        M[(iu[1], iu[0])] = sym_l

        # Diagonal = 1.0 for valid, NaN for invalid
        for i in range(n):
            M[i, i] = 1.0 if valid[i] else np.nan

        # Write matrix CSV (names as header and index)
        mat_df = pd.DataFrame(M, index=names, columns=names)
        mat_df.to_csv(args.out_matrix, index=True)
        print(f"[rank 0] Wrote matrix: {args.out_matrix} (shape {M.shape})")

        # Optional pairs TSV (upper triangle only, thresholded)
        if args.out_pairs is not None:
            thresh = float(args.pairs_thresh)
            out_path = args.out_pairs
            count = 0
            with open(out_path, "w") as f:
                f.write("i\tj\tname_i\tname_j\ttanimoto\n")
                for i in range(n):
                    # j > i to avoid duplicates & diagonal
                    row = M[i, (i+1):]
                    if np.all(np.isnan(row)):
                        continue
                    # Find indices passing threshold
                    for k, s in enumerate(row, start=i+1):
                        if not math.isnan(float(s)) and s >= thresh:
                            f.write(f"{i}\t{k}\t{names[i]}\t{names[k]}\t{float(s):.6f}\n")
                            count += 1
                            if count % 500000 == 0:
                                print(f"[rank 0] pairs written: {count}", flush=True)
            print(f"[rank 0] Wrote pairs: {out_path} (>= {thresh}) with {count} rows")

if __name__ == "__main__":
    main()
