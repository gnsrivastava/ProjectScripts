#!/bin/python3

import argparse, math
from typing import Optional, Dict, Tuple, List
import pandas as pd
from scipy.optimize import linear_sum_assignment

COLS = ["qseqid","sseqid","pident","length","qlen","slen",
        "qstart","qend","sstart","send","evalue","bitscore"]

def read_outfmt6(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", header=None, names=COLS,
                     usecols=["qseqid","sseqid","pident","length","evalue","bitscore"])
    for c in ["pident","length","evalue","bitscore"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=["length","evalue","bitscore"]).reset_index(drop=True)
    return df

def collapse_max_bitscore(df: pd.DataFrame) -> pd.DataFrame:
    """Keep the row with max bitscore for each (qseqid, sseqid)."""
    if df.empty:
        return df
    idx = df.groupby(["qseqid","sseqid"])["bitscore"].idxmax()
    return df.loc[idx].reset_index(drop=True)

def to_dir_dict(df: pd.DataFrame, flip: bool,
                emax: Optional[float], Lmin: Optional[int]) -> Dict[Tuple[str,str], dict]:
    """
    Return dict: (A,B) -> best row dict incl. bitscore/pident/length/evalue.
    If flip=True, store as (sseqid,qseqid) so both files key as (A,B).
    """
    D: Dict[Tuple[str,str], dict] = {}
    for _, r in df.iterrows():
        if Lmin is not None and r.length < Lmin:
            continue
        if emax  is not None and r.evalue > emax:
            continue
        A, B = (r.sseqid, r.qseqid) if flip else (r.qseqid, r.sseqid)
        key = (A, B)
        prev = D.get(key)
        if prev is None or r.bitscore > prev["bitscore"]:
            D[key] = {
                "qseqid": A,
                "sseqid": B,
                "bitscore": float(r.bitscore),
                "pident": float(r.pident) if pd.notna(r.pident) else float("nan"),
                "length": float(r.length),
                "evalue": float(r.evalue),
            }
    return D

def _avg_two(a: Optional[float], b: Optional[float]) -> float:
    """Average two numbers if both are finite; otherwise return the finite one, else NaN."""
    a_ok = a is not None and not pd.isna(a)
    b_ok = b is not None and not pd.isna(b)
    if a_ok and b_ok:
        return 0.5*(a+b)
    if a_ok:
        return a
    if b_ok:
        return b
    return float("nan")

def merge_symmetric(a2b: dict, b2a: dict, how: str = "avg"):
    """
    Merge directional dicts into a symmetric record per (A,B).
    Adds pident_a2b, pident_b2a, and avg_pident to each record.
    """
    merged = {}
    keys = set(a2b) | set(b2a)
    for k in keys:
        v1, v2 = a2b.get(k), b2a.get(k)  # v1 from A->B file, v2 from B->A file (flipped)
        if v1 and v2:
            if how == "max":
                chosen = v1 if v1["bitscore"] >= v2["bitscore"] else v2
                score = max(v1["bitscore"], v2["bitscore"])
            elif how == "min":
                chosen = v1 if v1["bitscore"] <= v2["bitscore"] else v2
                score = min(v1["bitscore"], v2["bitscore"])
            else:  # avg
                chosen = v1 if v1["bitscore"] >= v2["bitscore"] else v2
                score = 0.5*(v1["bitscore"] + v2["bitscore"])
            rec = chosen.copy()
            rec["score"] = score
            # add directional pidents
            rec["pident_a2b"] = float(v1["pident"])
            rec["pident_b2a"] = float(v2["pident"])
            rec["avg_pident"] = _avg_two(rec["pident_a2b"], rec["pident_b2a"])
            merged[k] = rec
        else:
            # Only one side present
            only = (v1 or v2).copy()
            only["score"] = only["bitscore"]
            if v1:  # came from A->B
                only["pident_a2b"] = float(v1["pident"])
                only["pident_b2a"] = float("nan")
            else:   # came from B->A
                only["pident_a2b"] = float("nan")
                only["pident_b2a"] = float(v2["pident"])
            only["avg_pident"] = _avg_two(only["pident_a2b"], only["pident_b2a"])
            merged[k] = only
    return merged

def build_cost_matrix(merged: dict):
    """Return queries, targets, score_matrix, cost_matrix (cost = max_score - score)."""
    queries = sorted({k[0] for k in merged})
    targets = sorted({k[1] for k in merged})
    qi = {q:i for i,q in enumerate(queries)}
    ti = {t:i for i,t in enumerate(targets)}
    n, m = len(queries), len(targets)
    score = [[0.0 for _ in range(m)] for __ in range(n)]
    cell = {}
    max_score = 0.0
    for (A,B), rec in merged.items():
        i, j = qi[A], ti[B]
        s = rec["score"]
        if s > score[i][j]:
            score[i][j] = s
            cell[(i,j)] = rec
        if s > max_score:
            max_score = s
    cost = [[(max_score - score[i][j]) for j in range(m)] for i in range(n)]
    return queries, targets, score, cost, cell, max_score

def solve_hungarian(cost, transpose_if_needed=True):
    """
    SciPy’s linear_sum_assignment minimizes over rows. If rows>cols, transpose so it can assign.
    Returns (row_idx, col_idx, transposed_flag).
    """
    n, m = len(cost), (len(cost[0]) if cost else 0)
    transposed = False
    if transpose_if_needed and n > m:
        # transpose to shape (m, n)
        tcost = list(map(list, zip(*cost)))
        col_ind, row_ind = linear_sum_assignment(tcost)  # swapped
        transposed = True
        return row_ind, col_ind, transposed
    else:
        row_ind, col_ind = linear_sum_assignment(cost)
        return row_ind, col_ind, transposed

def main():
    ap = argparse.ArgumentParser(description="Symmetric BLAST merge on bitscore, Hungarian via SciPy (with directional pidents).")
    ap.add_argument("a2b")
    ap.add_argument("b2a")
    ap.add_argument("-o","--out", default="assignment_bitscore.csv")
    ap.add_argument("--summary", default="summary_bitscore.txt")
    ap.add_argument("--how", choices=["avg","max","min"], default="avg",
                    help="merge both directions by avg/max/min bitscore")
    ap.add_argument("--emax", type=float, default=1e-3, help="keep HSPs with E <= emax")
    ap.add_argument("--Lmin", type=int, default=30, help="keep HSPs with alignment length >= Lmin")
    args = ap.parse_args()

    A = collapse_max_bitscore(read_outfmt6(args.a2b))
    B = collapse_max_bitscore(read_outfmt6(args.b2a))
    dA = to_dir_dict(A, flip=False, emax=args.emax, Lmin=args.Lmin)   # A->B
    dB = to_dir_dict(B, flip=True,  emax=args.emax, Lmin=args.Lmin)   # B->A (flipped to A,B keys)
    merged = merge_symmetric(dA, dB, how=args.how)

    queries, targets, score, cost, cell, max_score = build_cost_matrix(merged)
    r, c, transposed = solve_hungarian(cost, transpose_if_needed=True)

    # Collect chosen real pairs
    chosen = []
    n, m = len(queries), len(targets)
    for i, j in zip(r, c):
        ii, jj = (j, i) if transposed else (i, j)
        if ii < n and jj < m and (ii, jj) in cell:
            rec = cell[(ii, jj)]
            chosen.append({
                "qseqid": rec["qseqid"],
                "sseqid": rec["sseqid"],
                "bitscore": rec["bitscore"],
                "score": rec["score"],
                "pident_a2b": rec.get("pident_a2b", float("nan")),
                "pident_b2a": rec.get("pident_b2a", float("nan")),
                "avg_pident": rec.get("avg_pident", float("nan")),
                "length": rec["length"],
                "evalue": rec["evalue"],
            })

    out_df = pd.DataFrame(chosen)
    out_df.to_csv(args.out, index=False)

    with open(args.summary, "w") as fh:
        fh.write(f"Unique queries: {len(queries)}\n")
        fh.write(f"Unique targets: {len(targets)}\n")
        fh.write(f"Assignments chosen: {len(out_df)}\n")
        fh.write(f"Max bitscore across matrix: {max_score:.3f}\n")
        fh.write(f"Merge mode: {args.how}\n")
        fh.write(f"Filters: emax <= {args.emax}, Lmin >= {args.Lmin}\n")
        fh.write(f"Cost = max_bitscore - score (minimized)\n")
        fh.write(f"Transposed for LSA: {transposed}\n")

    print(f"Wrote {args.out} and {args.summary}")

if __name__ == "__main__":
    main()
