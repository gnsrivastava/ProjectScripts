from Bio import SeqIO
import os
import sys
import argparse

def resolve_in_path(filename: str, in_dir: str) -> str:
    # If user passes an absolute path or already a path, use it as-is
    return filename if os.path.isabs(filename) or os.path.dirname(filename) \
           else os.path.join(in_dir, filename)

def parse_fasta(input_filename: str,
                in_dir: str = "./faa_mod",
                out_dir: str = "./faa_cleaned",
                exclude_phrase: str = "hypothetical protein") -> None:
    input_path = resolve_in_path(input_filename, in_dir)
    print(input_path)

    if not os.path.isfile(input_path):
        sys.exit(f"[ERROR] Input FASTA not found: {input_path}")

    # ensure output directory exists
    os.makedirs(out_dir, exist_ok=True)

    sequences = SeqIO.parse(input_path, "fasta")
    needle = exclude_phrase.lower()
    filtered = [seq for seq in sequences if needle not in seq.description.lower()]

    output_path = os.path.join(out_dir, os.path.basename(input_filename))
    print(output_path)
    with open(output_path, "w") as output:
        SeqIO.write(filtered, output, "fasta")

    print(f"[OK] Wrote {len(filtered)} sequences to {output_path}")

def main():
    ap = argparse.ArgumentParser(
        description="Filter out sequences whose description contains a phrase "
                    "(default: 'hypothetical protein')."
    )
    ap.add_argument("-f", "--file", required=True,
                    help="FASTA filename (found in --in-dir) or full path.")
    ap.add_argument("--in-dir", default="./faa_mod",
                    help="Input directory (default: ./faa_mod).")
    ap.add_argument("--out-dir", default="./faa_cleaned",
                    help="Output directory (default: ./faa_cleaned).")
    ap.add_argument("--exclude", default="hypothetical protein",
                    help="Case-insensitive phrase to exclude (default: 'hypothetical protein').")
    args = ap.parse_args()

    parse_fasta(args.file, args.in_dir, args.out_dir, args.exclude)

if __name__ == "__main__":
    main()
