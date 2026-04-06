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
MIN_SPAN="$1"
MAX_TIME="$2"
CHR="$3"
TREE_DIR="$4"   
TRUE_DIR="$5"  
NUM_REPS="$6"
shift 6

# read MAPS until "--", then read POPS after
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

# decode task id to which map + replicate
tid="${SLURM_ARRAY_TASK_ID}"
map_idx=$(( (tid - 1) / NUM_REPS ))
rep_id=$(( (tid - 1) % NUM_REPS + 1 ))
if (( map_idx < 0 || map_idx >= ${#MAPS[@]} )); then
    echo "ERROR: task ${tid} produced map_idx=${map_idx}, but n_maps=${#MAPS[@]}" >&2
    exit 1
fi
gmap="${MAPS[$map_idx]}"
seed="${rep_id}"

# for each pop, use its pop-specific ts to find ibd segments
log_msg "Detected true ibd specifications
MIN_SPAN=${MIN_SPAN}
MAX_TIME=${MAX_TIME}
TREE_DIR=${TREE_DIR}
TRUE_DIR=${TRUE_DIR}
CHR=${CHR}"

mkdir -p "$TRUE_DIR"
for pop in "${POPS[@]}"; do
    in_ts="${TREE_DIR}/${gmap}_${seed}_${pop}.ts"
    out_tsv="${TRUE_DIR}/${gmap}_${seed}_${pop}.tsv"
    # check ts file was created
    if [[ -s "${in_ts}.tsz" ]]; then
        log_msg "found compressed ts (${in_ts}.tsz); unzipping"
        tsunzip "${in_ts}.tsz"
    elif [[ -s "${in_ts}" ]]; then
        log_msg "found uncompressed ts (${in_ts})"
    else
        echo "ERROR: missing input tree sequence (${in_ts} or ${in_ts}.tsz)" >&2
        exit 1
    fi

    log_msg "generating true ibd segments in ${out_tsv}"
    python /home1/karatas/ibd_compare/job_scripts/get_true_ibd.py \
        "$in_ts" "$out_tsv" "$CHR" "$MIN_SPAN" "$MAX_TIME"

    log_msg "rezipping input tree sequence to ${in_ts}.tsz"
    tszip "$in_ts"
    rm -f "$in_ts"
done

log_msg "done with ${gmap} seed=${seed}"
