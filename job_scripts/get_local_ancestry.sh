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

: "${SLURM_ARRAY_TASK_ID:?ERROR: run as a Slurm array}"

tree_dir="$1"
anc_dir="$2"
global_anc_dir="$3"
chr="$4"
sample_size="$5"
admixture_time="$6"
census_time_offset="$7"
num_reps="$8"
shift 8

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

all_outputs_exist=true
global_out_tsv="${global_anc_dir}/${gmap}_${seed}_ADMIX.tsv"
admix_found=false
for pop in "${pops[@]}"; do
    if [[ "${pop}" == "ADMIX" ]]; then
        admix_found=true
    fi
    out_tsv="${anc_dir}/${gmap}_${seed}_${pop}.tsv"
    if [[ ! -s "${out_tsv}" ]]; then
        all_outputs_exist=false
        break
    fi
done

if [[ "${admix_found}" != "true" ]]; then
    echo "ERROR: ADMIX population missing from pops list" >&2
    exit 1
fi

if [[ "${all_outputs_exist}" == "true" && -s "${global_out_tsv}" ]]; then
    log_msg "all ancestry outputs exist for ${gmap} seed=${seed}; skipping"
    exit 0
fi

in_ts="${tree_dir}/${gmap}_${seed}_all.ts"
temp_root="${TMPDIR:-/tmp}"
temp_dir="$(mktemp -d \
    "${temp_root}/local_anc.${gmap}_${seed}_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.XXXXXX")"
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

mkdir -p "${anc_dir}" "${global_anc_dir}"
log_msg "local ancestry task
task_id=${tid}
gmap=${gmap}
seed=${seed}
chr=${chr}
sample_size=${sample_size}
canonical_ts=${in_ts}
work_ts=${temp_ts}
pops=${pops[*]}
global_out_tsv=${global_out_tsv}"

admix_start=""
admix_end=""
for idx in "${!pops[@]}"; do
    pop="${pops[$idx]}"
    start=$(( idx * sample_size ))
    end=$(( (idx + 1) * sample_size - 1 ))
    out_tsv="${anc_dir}/${gmap}_${seed}_${pop}.tsv"

    if [[ "${pop}" == "ADMIX" ]]; then
        admix_start="${start}"
        admix_end="${end}"
    fi

    if [[ -s "${out_tsv}" ]]; then
        log_msg "local ancestry exists (${out_tsv}); skipping ${pop}"
        continue
    fi

    python /home1/karatas/ibd_compare/job_scripts/get_local_ancestry.py \
        "${temp_ts}" \
        "${out_tsv}" \
        "local" \
        "${chr}" \
        "${pop}" \
        "${start}" \
        "${end}" \
        "${admixture_time}" \
        "${census_time_offset}"
done

if [[ -z "${admix_start}" || -z "${admix_end}" ]]; then
    echo "ERROR: failed to locate ADMIX sample bounds" >&2
    exit 1
fi

if [[ -s "${global_out_tsv}" ]]; then
    log_msg "global ancestry exists (${global_out_tsv}); skipping ADMIX global"
else
    python /home1/karatas/ibd_compare/job_scripts/get_local_ancestry.py \
        "${temp_ts}" \
        "${global_out_tsv}" \
        "global" \
        "${chr}" \
        "ADMIX" \
        "${admix_start}" \
        "${admix_end}" \
        "${admixture_time}" \
        "${census_time_offset}"
fi

log_msg "done with local and global ancestry for ${gmap} seed=${seed}"
