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

TREE_DIR="$1"
ANC_DIR="$2"
CHR="$3"
SAMPLE_SIZE="$4"
ADMIXTURE_TIME="$5"
CENSUS_TIME_OFFSET="$6"
NUM_REPS="$7"
shift 7

MAPS=()
while [[ $# -gt 0 && "$1" != "--" ]]; do
    MAPS+=( "$1" )
    shift
done
shift 1
POPS=( "$@" )

if [[ ${#MAPS[@]} -eq 0 ]]; then
    echo "ERROR: no MAPS provided" >&2
    exit 1
fi
if [[ ${#POPS[@]} -eq 0 ]]; then
    echo "ERROR: no POPS provided" >&2
    exit 1
fi

tid="${SLURM_ARRAY_TASK_ID}"
map_idx=$(( (tid - 1) / NUM_REPS ))
rep_id=$(( (tid - 1) % NUM_REPS + 1 ))
if (( map_idx < 0 || map_idx >= ${#MAPS[@]} )); then
    echo "ERROR: task ${tid} produced map_idx=${map_idx}, but n_maps=${#MAPS[@]}" >&2
    exit 1
fi
gmap="${MAPS[$map_idx]}"
seed="${rep_id}"

in_ts="${TREE_DIR}/${gmap}_${seed}_all.ts"
if [[ -s "${in_ts}.tsz" ]]; then
    log_msg "found compressed ts (${in_ts}.tsz); unzipping"
    tsunzip "${in_ts}.tsz"
elif [[ -s "${in_ts}" ]]; then
    log_msg "found uncompressed ts (${in_ts})"
else
    echo "ERROR: missing input tree sequence (${in_ts} or ${in_ts}.tsz)" >&2
    exit 1
fi

mkdir -p "${ANC_DIR}"
log_msg "local ancestry task
TASK_ID=${tid}
GMAP=${gmap}
SEED=${seed}
CHR=${CHR}
SAMPLE_SIZE=${SAMPLE_SIZE}
POPS=${POPS[*]}"

for idx in "${!POPS[@]}"; do
    pop="${POPS[$idx]}"
    start=$(( idx * SAMPLE_SIZE ))
    end=$(( (idx + 1) * SAMPLE_SIZE - 1 ))
    out_tsv="${ANC_DIR}/${gmap}_${seed}_${pop}.tsv"

    python /home1/karatas/ibd_compare/job_scripts/get_local_ancestry.py \
        "${in_ts}" \
        "${out_tsv}" \
        "${CHR}" \
        "${pop}" \
        "${start}" \
        "${end}" \
        "${ADMIXTURE_TIME}" \
        "${CENSUS_TIME_OFFSET}"
done

log_msg "rezipping input tree sequence to ${in_ts}.tsz"
tszip "${in_ts}"
rm -f "${in_ts}"

log_msg "done with local ancestry for ${gmap} seed=${seed}"
