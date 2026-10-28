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
#SBATCH --output=/home1/karatas/logs/trueIBD/%A_%a.%x.out
#SBATCH --error=/home1/karatas/logs/trueIBD/%A_%a.%x.err
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
ooa_tree_dir="$4"
true_dir="$5"  
sample_size="$6"
analysis_sample_size="$7"
num_reps="$8"
contract_file="$9"
id_contract="${10}"
shift 10
source /home1/karatas/ibd_compare/other_scripts/require_sim_id_contract.sh \
    "$contract_file" "$id_contract"

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
rep="${rep_id}"

log_msg "Detected true ibd specifications
min_span=${min_span}
max_time=${max_time}
ooa_tree_dir=${ooa_tree_dir}
true_dir=${true_dir}
chr=${chr}"

mkdir -p "$true_dir"
in_ts="${ooa_tree_dir}/${gmap}_${rep}_chr${chr}_all.ts"
temp_root="${TMPDIR:-/tmp}"
temp_dir="$(mktemp -d \
    "${temp_root}/true_ibd.${gmap}_${rep}_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.XXXXXX")"
temp_ts="${temp_dir}/${gmap}_${rep}_chr${chr}_all.ts"
cleanup() {
    rm -rf "${temp_dir}"
}
trap cleanup EXIT

if [[ -s "${in_ts}.tsz" ]]; then
    temp_tsz="${temp_dir}/${gmap}_${rep}_chr${chr}_all.ts.tsz"
    cp "${in_ts}.tsz" "${temp_tsz}"
    log_msg "staging compressed ts
canonical_input=${in_ts}.tsz
temp_input=${temp_tsz}
temp_output=${temp_ts}"
    tsunzip "${temp_tsz}"
else
    echo "ERROR: missing OOA tree sequence ${in_ts}.tsz" >&2
    exit 1
fi

for idx in "${!pops[@]}"; do
    pop="${pops[$idx]}"
    start=$(( idx * sample_size ))
    end=$(( start + analysis_sample_size - 1 ))
    out_tsv="${true_dir}/${gmap}_${rep}_chr${chr}_${pop}.tsv"

    log_msg "generating true ibd segments in ${out_tsv}"
    python /home1/karatas/ibd_compare/job_scripts/get_true_ibd.py \
        "$temp_ts" "$out_tsv" "$chr" "$min_span" "$max_time" "$start" "$end" \
        "$pop"
done

log_msg "done with ${gmap} rep=${rep}"
