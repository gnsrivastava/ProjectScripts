#!/bin/bash
#Note: Soft masking helps produce more biologically meaningful results by reducing the number of spurious, or false-positive, matches that arise from common, non-specific sequence patterns.

#BLASTP=${BLASTP:-blastp} # I have the location of the blastp in source file. Please uncomment according to your needs. 
THREADS_PER_JOB=${THREADS_PER_JOB:-1}
JOBS=${JOBS:-20}
OUTDIR=${OUTDIR:-blast_out1}
SPECIES_FILE=${SPECIES_FILE:-species.txt} # list of the species with one species per row. Used to create species pairs

mkdir -p "$OUTDIR"

run_one() {
  local A="$1" B="$2"
  local q="./MIC_Analysis/${A}_nonredundant.faa"
  local db="./MIC_Analysis/${B}/${B}_db"
  local out="${OUTDIR}/${A}_vs_${B}.tsv"

  if [[ -s "$out" ]]; then
    echo "SKIP: ${A}_vs_${B}"
    return 0
  fi

  echo "RUN : ${A}_vs_${B}"
  "$BLASTP" \
    -db "$db" \
    -query "$q" \
    -task blastp \
    -matrix BLOSUM62 \
    -comp_based_stats 2 \
    -seg yes -soft_masking true \ # a filtering method used to disregard low-complexity regions during the initial, computationally intensive stages of a BLAST search
    -gapopen 11 -gapextend 1 \
    -max_hsps 1 \
    -evalue 1e6 \
    -use_sw_tback \
    -outfmt '6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore' \
  > "$out"
}
export -f run_one
export BLASTP THREADS_PER_JOB OUTDIR

# ---- Read species file without process substitution ----
if [[ ! -f "$SPECIES_FILE" ]]; then
  echo "ERROR: ${SPECIES_FILE} not found" >&2
  exit 1
fi

# Build an array SPECIES (skip blank/whitespace-only lines)
SPECIES=()
while IFS= read -r line; do
  [[ -z "${line//[[:space:]]/}" ]] && continue
  SPECIES+=("$line")
done < "$SPECIES_FILE"

# Generate ordered pairs A!=B
pairs_file="$(mktemp)"
trap 'rm -f "$pairs_file"' EXIT
for A in "${SPECIES[@]}"; do
  for B in "${SPECIES[@]}"; do
    [[ "$A" == "$B" ]] && continue
    printf "%s\t%s\n" "$A" "$B" >> "$pairs_file"
  done
done

# Run
if command -v parallel >/dev/null 2>&1; then
  parallel --colsep '\t' -j "$JOBS" run_one {1} {2} :::: "$pairs_file"
else
  while IFS=$'\t' read -r A B; do
    run_one "$A" "$B"
  done < "$pairs_file"
fi

echo "All done. Results in: $OUTDIR"
