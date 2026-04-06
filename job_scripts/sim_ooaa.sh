#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           sim_nonpyrho.sh
###############################################################################

#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --partition=qcb
#SBATCH --account=jazlynmo_738
#SBATCH --nodes=1
#SBATCH --output=/home1/karatas/logs/sim.%A_%a.%x.out
#SBATCH --error=/home1/karatas/logs/sim.%A_%a.%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=karatas@usc.edu

set -euo pipefail
module purge
ml gcc/13.3.0 htslib/1.19.1 conda
source /apps/conda/miniforge3/25.3.0/etc/profile.d/conda.sh
conda activate stdpopsim
export PATH="/home1/karatas/.conda/envs/stdpopsim/bin:$PATH"
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh

# constants
: "${SLURM_ARRAY_TASK_ID:?ERROR: run as a Slurm array}"

# inputs
TREE_DIR="$1"
SAMPLE_SIZE="$2"
NUM_REPS="$3"
MSPRIME_MODEL="$4"
DEMO_MODEL="$5"
CHR="$6"
shift 6
# decode task id to which map + replicate
MAPS=( "$@" )
if [[ ${#MAPS[@]} -eq 0 ]]; then
    echo "ERROR: no GMAPs provided (after shifting args)" >&2
    exit 1
fi
tid="${SLURM_ARRAY_TASK_ID}"
map_idx=$(( (tid - 1) / NUM_REPS ))      # 0..n_maps-1
rep_id=$(( (tid - 1) % NUM_REPS + 1 ))   # 1..NUM_REPS
if (( map_idx < 0 || map_idx >= ${#MAPS[@]} )); then
    echo "ERROR: task ${tid} produced map_idx=${map_idx}, but n_maps=${#MAPS[@]}" >&2
    exit 1
fi
gmap="${MAPS[$map_idx]}"
seed="${rep_id}"

mkdir -p "${TREE_DIR}"

log_msg "stdpopsim sim (combined array)
TASK_ID=${tid}
GMAP=${gmap} (map_idx=${map_idx}/${#MAPS[@]})
SEED=${seed}
TREE_DIR=${TREE_DIR}
SAMPLE_SIZE=${SAMPLE_SIZE}
MODEL=${MSPRIME_MODEL}
sample_spec=${sample_spec[*]}"

# TODO:
python /home1/karatas/ibd_compare/job_scripts/python_util/sim_ooaa.py \
        
# compress tree sequence
log_msg "compress tree sequence"
tszip "$out_ts"
rm -f "$out_ts"

log_msg "sim done"
