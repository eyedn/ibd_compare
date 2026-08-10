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
#SBATCH --output=/home1/karatas/logs/trueFF/%A_%a.%x.out
#SBATCH --error=/home1/karatas/logs/trueFF/%A_%a.%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=karatas@usc.edu

set -euo pipefail
module purge
ml rstats/4.5.1
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh

: "${SLURM_ARRAY_TASK_ID:?ERROR: run as a Slurm array}"

# inputs
true_dir="$1"
common_ff_dir="$2"
ref_dir="$3"
chr="$4"
min_bp_len="$5" 
min_cm_len="$6" 
num_reps="$7"
contract_file="$8"
id_contract="$9"
shift 9
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

mkdir -p "${common_ff_dir}/true_ibd"

# loop through each pop
for pop in "${pops[@]}"; do
    true_ibd="${true_dir}/${gmap}_${rep}_chr${chr}_${pop}.tsv"

    # determine gmap path for interpolation; create standard gmap labels
    if [[ "$gmap" == *HapMapII_GRCh38* ]]; then 
        gmap_file="${ref_dir}/HapMapII_GRCh38/plink.chr${chr}.GRCh38.UCSC.format.map"
        gmap_label="hapmap"
    elif [[ "$gmap" == *DeCodeSexAveraged_GRCh38* ]]; then
        gmap_file="${ref_dir}/deCODE/addendum/sex-averaged_gmaps/sex_av_chr${chr}_hg38.gmap"
        gmap_label="decode"
    elif [[ "$gmap" == *DefaultUniform* ]]; then
        gmap_file="${ref_dir}/GRCh38.p14/recomb_maps_flat1.15235e-08/chr${chr}.map"
        gmap_label="flat"
    fi
    
    # run formatting
    # note, "-1" set for inconsequential variables for formatting true ibd
    ff_out="${common_ff_dir}/true_ibd/true_${gmap_label}_${gmap_label}_${pop}_${rep}.bed"
    log_msg "formatting true segments for ${gmap_label} pop${pop} rep${rep}"
    Rscript "/home1/karatas/ibd_compare/job_scripts/filter_and_format.R" \
        "$min_bp_len" \
        "$min_cm_len" \
        "-1" \
        "-1" \
        "$true_ibd" \
        "$gmap_file" \
        "$ff_out" \
        "-1" \
        "TRUE" \
        "FALSE"
    gzip -f "$ff_out"
done
