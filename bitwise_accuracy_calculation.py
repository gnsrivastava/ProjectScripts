import pandas as pd
import numpy as np

# -------------------------
# Helpers
# -------------------------
def parse_species(protein_id: str) -> str:
    # "1006581.GCW_00035" -> "1006581"
    return str(protein_id).split(".", 1)[0]

def split_ec_list(ec_str: str):
    # "1.1.5.3,1.4.3.19" -> ["1.1.5.3","1.4.3.19"]
    if pd.isna(ec_str):
        return []
    s = str(ec_str).strip()
    if not s or s.lower() == "none":
        return []
    # allow commas or semicolons
    parts = [p.strip() for p in s.replace(";", ",").split(",")]
    return [p for p in parts if p]

def clean_ec(ec: str) -> str:
    # remove "EC:" and spaces; keep '-' if present
    return str(ec).replace("EC:", "").replace("EC ", "").strip()

def ec_levels(ec: str):
    # "3.1.26.5" -> ["3","1","26","5"]
    # if incomplete or malformed, return []
    ec = clean_ec(ec)
    if "-" in ec:  # treat ambiguous as prefix match only; levels from prefix part
        ec = ec.split("-", 1)[0]
    parts = ec.split(".")
    if len(parts) < 1:
        return []
    return parts

def bitwise_counts_for_pair(pred_ec: str, true_ec: str):
    """
    Returns (c1,c2,c3,c4) for this pred/true pair:
    counts how many prefix levels match consecutively from level 1.
    """
    pred = ec_levels(pred_ec)
    true = ec_levels(true_ec)
    if not pred or not true:
        return (0, 0, 0, 0)

    maxk = min(4, len(pred), len(true))
    k = 0
    for i in range(maxk):
        if pred[i] == true[i]:
            k += 1
        else:
            break

    return (1, 0, 0, 0) if k == 1 else \
           (1, 1, 0, 0) if k == 2 else \
           (1, 1, 1, 0) if k == 3 else \
           (1, 1, 1, 1) if k >= 4 else \
           (0, 0, 0, 0)

def best_bitwise_counts(pred_ecs, true_ecs):
    """
    For one protein that can have multiple true ECs and one predicted EC (or multiple),
    return the best matching prefix counts across all pred×true combinations.
    """
    best = (0, 0, 0, 0)
    for p in pred_ecs:
        for t in true_ecs:
            c = bitwise_counts_for_pair(p, t)
            # choose the combination with highest number of matched levels
            if sum(c) > sum(best):
                best = c
            if best == (1, 1, 1, 1):
                return best
    return best

# -------------------------
# Main: per-species bitwise accuracy
# -------------------------
def bitwise_accuracy_by_species(merged_df: pd.DataFrame) -> pd.DataFrame:
    """
    Expects columns:
      - protein
      - EC_number (true; can be comma-separated)
      - EC_4d_label (predicted full EC as string; can be comma-separated, but usually one)
    Computes per-species mean bitwise accuracies (Bit1..Bit4) across proteins.
    """
    df = merged_df.copy()

    # Species from protein
    df["species"] = df["protein"].map(parse_species)

    # True ECs list
    df["true_ecs"] = df["EC_number"].map(split_ec_list)

    # Predicted ECs list (use EC_4d_label as the predicted EC)
    df["pred_ecs"] = df["EC_4d_label"].map(split_ec_list)

    # Per-protein counts
    c1 = np.zeros(len(df), dtype=np.int64)
    c2 = np.zeros(len(df), dtype=np.int64)
    c3 = np.zeros(len(df), dtype=np.int64)
    c4 = np.zeros(len(df), dtype=np.int64)

    for i, (preds, trues) in enumerate(zip(df["pred_ecs"], df["true_ecs"])):
        b1, b2, b3, b4 = best_bitwise_counts(preds, trues)
        c1[i], c2[i], c3[i], c4[i] = b1, b2, b3, b4

    df["bit1"] = c1
    df["bit2"] = c2
    df["bit3"] = c3
    df["bit4"] = c4

    # Per-species accuracy = mean of bits across proteins in that species
    out = (
        df.groupby("species", as_index=False)
          .agg(
              n_proteins=("protein", "count"),
              Bit1=("bit1", "mean"),
              Bit2=("bit2", "mean"),
              Bit3=("bit3", "mean"),
              Bit4=("bit4", "mean"),
          )
          .sort_values("species")
    )
    return out

# -------------------------
# Example usage
# -------------------------
merged = pd.read_csv('merged_data_protein.tsv', sep='\t')
species_acc = bitwise_accuracy_by_species(merged)
species_acc.to_csv("bitwise_accuracy_by_species.tsv", sep="\t", index=False)
