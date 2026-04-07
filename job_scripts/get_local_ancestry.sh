#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           get_local_ancestry.sh
###############################################################################


set -euo pipefail

module purge
ml conda
source /apps/conda/miniforge3/25.3.0/etc/profile.d/conda.sh
conda activate stdpopsim
export PATH="/home1/karatas/.conda/envs/stdpopsim/bin:$PATH"
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh

: "${SLURM_ARRAY_TASK_ID:?ERROR: run as a Slurm array}"

tree_dir="$1"
anc_dir="$2"
chr="$3"
sample_size="$4"
admixture_time="$5"
census_time_offset="$6"
num_reps="$7"
shift 7

maps=()
while [[ $# -gt 0 && "$1" != "--" ]]; do
    maps+=( "$1" )
    shift
done
shift 1
pops=( "$@" )

if [[ ${#maps[@]} -eq 0 ]]; then
    echo "ERROR: no maps provided" >&2
    exit 1
fi
if [[ ${#pops[@]} -eq 0 ]]; then
    echo "ERROR: no pops provided" >&2
    exit 1
fi

tid="${SLURM_ARRAY_TASK_ID}"
map_idx=$(( (tid - 1) / num_reps ))
rep_id=$(( (tid - 1) % num_reps + 1 ))
if (( map_idx < 0 || map_idx >= ${#maps[@]} )); then
    echo "ERROR: task ${tid} produced map_idx=${map_idx}, but n_maps=${#maps[@]}" >&2
    exit 1
fi
gmap="${maps[$map_idx]}"
seed="${rep_id}"

in_ts="${tree_dir}/${gmap}_${seed}_all.ts"
if [[ -s "${in_ts}.tsz" ]]; then
    log_msg "found compressed ts (${in_ts}.tsz); unzipping"
    tsunzip "${in_ts}.tsz"
elif [[ -s "${in_ts}" ]]; then
    log_msg "found uncompressed ts (${in_ts})"
else
    echo "ERROR: missing input tree sequence (${in_ts} or ${in_ts}.tsz)" >&2
    exit 1
fi

mkdir -p "${anc_dir}"
log_msg "local ancestry task
task_id=${tid}
gmap=${gmap}
seed=${seed}
chr=${chr}
sample_size=${sample_size}
pops=${pops[*]}"

for idx in "${!pops[@]}"; do
    pop="${pops[$idx]}"
    start=$(( idx * sample_size ))
    end=$(( (idx + 1) * sample_size - 1 ))
    out_tsv="${anc_dir}/${gmap}_${seed}_${pop}.tsv"

    python /home1/karatas/ibd_compare/job_scripts/get_local_ancestry.py \
        "${in_ts}" \
        "${out_tsv}" \
        "${chr}" \
        "${pop}" \
        "${start}" \
        "${end}" \
        "${admixture_time}" \
        "${census_time_offset}"
done

log_msg "rezipping input tree sequence to ${in_ts}.tsz"
tszip "${in_ts}"
rm -f "${in_ts}"

log_msg "done with local ancestry for ${gmap} seed=${seed}"
