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
#SBATCH --output=/home1/karatas/logs/combTVI/%A_%a.%x.out
#SBATCH --error=/home1/karatas/logs/combTVI/%A_%a.%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=karatas@usc.edu

set -euo pipefail
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh

# input
tvi_dir="$1"
out="$2"
out_per_rep="$3"
contract_file="$4"
id_contract="$5"

source /home1/karatas/ibd_compare/other_scripts/require_sim_id_contract.sh \
    "$contract_file" "$id_contract"

mkdir -p "$(dirname "$out")" "$(dirname "$out_per_rep")"

# combine aggregate summaries
files=()
while IFS= read -r file; do
    files+=( "$file" )
done < <(
    find "$tvi_dir" -type f -name "*.agg.tsv" ! -path "$out" \
        ! -path "$out_per_rep" | LC_ALL=C sort
)
[[ ${#files[@]} -gt 0 ]] || { log_msg "ERROR: No .agg.tsv files found in $tvi_dir" >&2; exit 1; }

tmp_out="$(mktemp "${out}.tmp.XXXXXX")"
trap 'rm -f "$tmp_out" "${tmp_per:-}"' EXIT
head -n 1 "${files[0]}" > "$tmp_out"

for f in "${files[@]}"; do
    cmp -s <(head -n 1 "${files[0]}") <(head -n 1 "$f") || {
        log_msg "ERROR: aggregate schema mismatch: $f" >&2
        exit 1
    }
    tail -n +2 "$f" >> "$tmp_out"
done

mv -f "$tmp_out" "$out"
log_msg "Wrote combined results: $out"
log_msg "Combined N(agg)=$((${#files[@]})) files"

# combine per-replicate summaries
per_rep_files=()
while IFS= read -r file; do
    per_rep_files+=( "$file" )
done < <(
    find "$tvi_dir" -type f -name "*.per_rep.tsv" ! -path "$out" \
        ! -path "$out_per_rep" | LC_ALL=C sort
)
[[ ${#per_rep_files[@]} -gt 0 ]] || { log_msg "ERROR: No *.per_rep.tsv files found in $tvi_dir" >&2; exit 1; }

tmp_per="$(mktemp "${out_per_rep}.tmp.XXXXXX")"
head -n 1 "${per_rep_files[0]}" > "$tmp_per"

for f in "${per_rep_files[@]}"; do
    cmp -s <(head -n 1 "${per_rep_files[0]}") <(head -n 1 "$f") || {
        log_msg "ERROR: per-replicate schema mismatch: $f" >&2
        exit 1
    }
    tail -n +2 "$f" >> "$tmp_per"
done

mv -f "$tmp_per" "$out_per_rep"
log_msg "Wrote combined per-replicate results: $out_per_rep"
log_msg "Combined N(per_rep)=$((${#per_rep_files[@]})) files"
trap - EXIT
