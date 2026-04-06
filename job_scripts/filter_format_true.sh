#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           filter_format_true.sh
###############################################################################

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --time=1-00:00:00
#SBATCH --partition=qcb
#SBATCH --account=jazlynmo_738
#SBATCH --nodes=1
#SBATCH --output=/home1/karatas/logs/trueFF.%A_%a.out
#SBATCH --error=/home1/karatas/logs/trueFF.%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=karatas@usc.edu

set -euo pipefail
module purge
ml rstats/4.5.1
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh

: "${SLURM_ARRAY_TASK_ID:?ERROR: run as a Slurm array}"

# inputs
TRUE_DIR="$1"
COMMON_FF_DIR="$2"
REF_DIR="$3"
CHR="$4"
MIN_BP_LEN="$5" 
MIN_CM_LEN="$6" 
NUM_REPS="$7"
shift 7

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

mkdir -p "${COMMON_FF_DIR}/true_ibd"

# loop through each pop
for pop in "${POPS[@]}"; do
    true_ibd="${TRUE_DIR}/${gmap}_${seed}_${pop}.tsv"

    # determine gmap path for interpolation; create standard gmap labels
    if [[ "$gmap" == *HapMapII_GRCh38* ]]; then 
        gmap_file="${REF_DIR}/HapMapII_GRCh38/plink.chr${CHR}.GRCh38.UCSC.format.map"
        gmap_label="hapmap"
    elif [[ "$gmap" == *DeCodeSexAveraged_GRCh38* ]]; then
        gmap_file="${REF_DIR}/deCODE/addendum/sex-averaged_gmaps/sex_av_chr${CHR}_hg38.gmap"
        gmap_label="decode"
    elif [[ "$gmap" == *DefaultUniform* ]]; then
        gmap_file="${REF_DIR}/GRCh38.p14/recomb_maps_flat1.15235e-08/chr${CHR}.map"
        gmap_label="flat"
    fi
    
    # run formatting
    # note, "-1" set for inconsequential variables for formatting true ibd
    ff_out="${COMMON_FF_DIR}/true_ibd/true_${gmap_label}_${gmap_label}_${pop}_${seed}.bed"
    log_msg "formatting true segments for ${gmap_label} pop${pop} rep${seed}"
    Rscript "/home1/karatas/ibd_compare/job_scripts/filter_and_format.R" \
        "$MIN_BP_LEN" \
        "$MIN_CM_LEN" \
        "-1" \
        "-1" \
        "$true_ibd" \
        "$gmap_file" \
        "$ff_out" \
        "-1" \
        "TRUE"
    gzip -f "$ff_out"
done
