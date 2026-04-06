#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           combine_tvi.sh
###############################################################################

#SBATCH -p qcb
#SBATCH -A jazlynmo_738
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH --output=/home1/karatas/logs/combTVI.%A_%a.%x.out
#SBATCH --error=/home1/karatas/logs/combTVI.%A_%a.%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=karatas@usc.edu

set -euo pipefail
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh

# input
TVI_DIR="$1"
OUT="$2"
OUT_PER_REP="$3"

mkdir -p "$(dirname "$OUT")" "$(dirname "$OUT_PER_REP")"

# combine aggregate summaries
mapfile -t files < <(
    find "$TVI_DIR" -type f -name "*.agg.tsv" | LC_ALL=C sort
)
[[ ${#files[@]} -gt 0 ]] || { log_msg "ERROR: No .agg.tsv files found in $TVI_DIR" >&2; exit 1; }

tmp_out="$(mktemp "${OUT}.tmp.XXXXXX")"
head -n 1 "${files[0]}" > "$tmp_out"

for f in "${files[@]}"; do
    tail -n +2 "$f" >> "$tmp_out"
done

mv -f "$tmp_out" "$OUT"
log_msg "Wrote combined results: $OUT"
log_msg "Combined N(agg)=$((${#files[@]})) files"

# combine per-replicate summaries
mapfile -t per_rep_files < <(
    find "$TVI_DIR" -type f -name "*.per_rep.tsv" | LC_ALL=C sort
)
[[ ${#per_rep_files[@]} -gt 0 ]] || { log_msg "ERROR: No *.per_rep.tsv files found in $TVI_DIR" >&2; exit 1; }

tmp_per="$(mktemp "${OUT_PER_REP}.tmp.XXXXXX")"
head -n 1 "${per_rep_files[0]}" > "$tmp_per"

for f in "${per_rep_files[@]}"; do
    tail -n +2 "$f" >> "$tmp_per"
done

mv -f "$tmp_per" "$OUT_PER_REP"
log_msg "Wrote combined per-replicate results: $OUT_PER_REP"
log_msg "Combined N(per_rep)=$((${#per_rep_files[@]})) files"
