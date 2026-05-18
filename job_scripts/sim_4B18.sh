#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           sim_4B18.sh
###############################################################################

# SBATCH commands if appropriate that were not specified in another
# script that directly calls this file

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
tree_dir="$1"
sample_size="$2"
num_reps="$3"
msprime_model="$4"
stdpopsim_demography_id="$5"
chr="$6"
shift 6
# decode task id to which map + replicate
maps=( "$@" )
if [[ ${#maps[@]} -eq 0 ]]; then
    echo "ERROR: no GMAPs provided (after shifting args)" >&2
    exit 1
fi
tid="${SLURM_ARRAY_TASK_ID}"
map_idx=$(( (tid - 1) / num_reps ))      # 0..n_maps-1
rep_id=$(( (tid - 1) % num_reps + 1 ))   # 1..num_reps
if (( map_idx < 0 || map_idx >= ${#maps[@]} )); then
    echo "ERROR: task ${tid} produced map_idx=${map_idx}, but n_maps=${#maps[@]}" >&2
    exit 1
fi
gmap="${maps[$map_idx]}"
seed="${rep_id}"

# check if this task was already completed; if so, skip
out_ts="${tree_dir}/${gmap}_${seed}_all.ts"
if [[ -s "${out_ts}.tsz" ]]; then
    log_msg "outputs exist (${out_ts}.tsz); skipping"
    exit 0
fi

# define sampling specification
sample_spec=(
    "AFR:${sample_size}"
    "EUR:${sample_size}"
    "ASIA:${sample_size}"
    "ADMIX:${sample_size}"
)

mkdir -p "${tree_dir}"

log_msg "4B18 stdpopsim sim (combined array)
task_id=${tid}
gmap=${gmap} (map_idx=${map_idx}/${#maps[@]})
seed=${seed}
tree_dir=${tree_dir}
sample_size=${sample_size}
model=${msprime_model}
demography_id=${stdpopsim_demography_id}
sample_spec=${sample_spec[*]}"

# run stdpopsim (adjust specifications here)
if [[ "${gmap}" == *DefaultUniform* ]]; then
    stdpopsim -v -e msprime --msprime-model "$msprime_model" \
        HomSap -s "$seed" -c "chr${chr}" \
        -o "$out_ts" \
        -d "$stdpopsim_demography_id" \
        "${sample_spec[@]}"
else
    stdpopsim -v -e msprime --msprime-model "$msprime_model" \
        HomSap -s "$seed" -g "$gmap" -c "chr${chr}" \
        -o "$out_ts" \
        -d "$stdpopsim_demography_id" \
        "${sample_spec[@]}"
fi

# compress tree sequence
log_msg "compress tree sequence"
tszip "$out_ts"
rm -f "$out_ts"

log_msg "sim done"
