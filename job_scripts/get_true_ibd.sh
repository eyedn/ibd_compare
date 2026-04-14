#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           get_true_ibd.sh
###############################################################################

#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=1-00:00:00
#SBATCH --partition=qcb
#SBATCH --account=jazlynmo_738
#SBATCH --nodes=1
#SBATCH --output=/home1/karatas/logs/slurm.%A_%a.%x.out
#SBATCH --error=/home1/karatas/logs/slurm.%A_%a.%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=karatas@usc.edu

set -euo pipefail

module purge
ml conda
source /apps/conda/miniforge3/25.3.0/etc/profile.d/conda.sh
conda activate stdpopsim
export PATH="/home1/karatas/.conda/envs/stdpopsim/bin:$PATH"
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh
log_msg "using python: $(which python)"

: "${SLURM_ARRAY_TASK_ID:?ERROR: run as a Slurm array}"

# inputs
min_span="$1"
max_time="$2"
chr="$3"
tree_dir="$4"   
true_dir="$5"  
sample_size="$6"
num_reps="$7"
shift 7

# read maps until "--", then read pops after
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

# decode task id to which map + replicate
tid="${SLURM_ARRAY_TASK_ID}"
map_idx=$(( (tid - 1) / num_reps ))
rep_id=$(( (tid - 1) % num_reps + 1 ))
if (( map_idx < 0 || map_idx >= ${#maps[@]} )); then
    echo "ERROR: task ${tid} produced map_idx=${map_idx}, but n_maps=${#maps[@]}" >&2
    exit 1
fi
gmap="${maps[$map_idx]}"
seed="${rep_id}"

log_msg "Detected true ibd specifications
min_span=${min_span}
max_time=${max_time}
tree_dir=${tree_dir}
true_dir=${true_dir}
chr=${chr}"

mkdir -p "$true_dir"
in_ts="${tree_dir}/${gmap}_${seed}_all.ts"
temp_root="${TMPDIR:-/tmp}"
temp_dir="$(mktemp -d \
    "${temp_root}/true_ibd.${gmap}_${seed}_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.XXXXXX")"
temp_ts="${temp_dir}/${gmap}_${seed}_all.ts"
cleanup() {
    rm -rf "${temp_dir}"
}
trap cleanup EXIT

if [[ -s "${in_ts}.tsz" ]]; then
    temp_tsz="${temp_dir}/${gmap}_${seed}_all.ts.tsz"
    cp "${in_ts}.tsz" "${temp_tsz}"
    log_msg "staging compressed ts
canonical_input=${in_ts}.tsz
temp_input=${temp_tsz}
temp_output=${temp_ts}"
    tsunzip "${temp_tsz}"
elif [[ -s "${in_ts}" ]]; then
    log_msg "staging uncompressed ts
canonical_input=${in_ts}
temp_output=${temp_ts}"
    cp "${in_ts}" "${temp_ts}"
else
    echo "ERROR: missing input tree sequence (${in_ts} or ${in_ts}.tsz)" >&2
    exit 1
fi

for idx in "${!pops[@]}"; do
    pop="${pops[$idx]}"
    start=$(( idx * sample_size ))
    end=$(( (idx + 1) * sample_size - 1 ))
    out_tsv="${true_dir}/${gmap}_${seed}_${pop}.tsv"

    log_msg "generating true ibd segments in ${out_tsv}"
    python /home1/karatas/ibd_compare/job_scripts/get_true_ibd.py \
        "$temp_ts" "$out_tsv" "$chr" "$min_span" "$max_time" "$start" "$end"
done

log_msg "done with ${gmap} seed=${seed}"
