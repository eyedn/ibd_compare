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


set -euo pipefail
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh

# input
tvi_dir="$1"
out="$2"
out_per_rep="$3"

mkdir -p "$(dirname "$out")" "$(dirname "$out_per_rep")"

# combine aggregate summaries
mapfile -t files < <(
    find "$tvi_dir" -type f -name "*.agg.tsv" | LC_ALL=C sort
)
[[ ${#files[@]} -gt 0 ]] || { log_msg "ERROR: No .agg.tsv files found in $tvi_dir" >&2; exit 1; }

tmp_out="$(mktemp "${out}.tmp.XXXXXX")"
head -n 1 "${files[0]}" > "$tmp_out"

for f in "${files[@]}"; do
    tail -n +2 "$f" >> "$tmp_out"
done

mv -f "$tmp_out" "$out"
log_msg "Wrote combined results: $out"
log_msg "Combined N(agg)=$((${#files[@]})) files"

# combine per-replicate summaries
mapfile -t per_rep_files < <(
    find "$tvi_dir" -type f -name "*.per_rep.tsv" | LC_ALL=C sort
)
[[ ${#per_rep_files[@]} -gt 0 ]] || { log_msg "ERROR: No *.per_rep.tsv files found in $tvi_dir" >&2; exit 1; }

tmp_per="$(mktemp "${out_per_rep}.tmp.XXXXXX")"
head -n 1 "${per_rep_files[0]}" > "$tmp_per"

for f in "${per_rep_files[@]}"; do
    tail -n +2 "$f" >> "$tmp_per"
done

mv -f "$tmp_per" "$out_per_rep"
log_msg "Wrote combined per-replicate results: $out_per_rep"
log_msg "Combined N(per_rep)=$((${#per_rep_files[@]})) files"
