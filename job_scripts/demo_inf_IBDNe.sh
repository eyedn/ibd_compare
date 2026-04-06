#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           demo_inf_GONE2.sh
###############################################################################

#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --partition=qcb
#SBATCH --account=jazlynmo_738
#SBATCH --nodes=1
#SBATCH --output=/home1/karatas/logs/ibdne.%A_%a.%x.out
#SBATCH --error=/home1/karatas/logs/ibdne.%A_%a.%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=karatas@usc.edu

set -euo pipefail
ml purge
ml gcc/13.3.0 htslib/1.19.1 bcftools/1.19

: "${SLURM_ARRAY_TASK_ID:?ERROR: run as a Slurm array}"
THREADS="${SLURM_CPUS_PER_TASK:-1}"

# inputs
COMMON_FF_DIR="$1"
DEMO_DIR="$2"
CHR="$3"
FLAT_RATE="$4"
NUM_REPS="$5"
shift 5

# read TMAPS until "--", then read POPS after
TMAPS=()
while [[ $# -gt 0 && "$1" != "--" ]]; do
    TMAPS+=( "$1" )
    shift
done
if [[ $# -eq 0 ]]; then
    echo "ERROR: expected '--' separator between TMAPS and POPS" >&2
    exit 1
fi
shift 1
POPS=( "$@" )

if [[ ${#TMAPS[@]} -eq 0 ]]; then
    echo "ERROR: no TMAPS provided" >&2
    exit 1
fi
if [[ ${#POPS[@]} -eq 0 ]]; then
    echo "ERROR: no POPS provided" >&2
    exit 1
fi

tid="${SLURM_ARRAY_TASK_ID}"
tmap_idx=$(( (tid - 1) / NUM_REPS ))
rep_id=$(( (tid - 1) % NUM_REPS + 1 ))

if (( tmap_idx < 0 || tmap_idx >= ${#TMAPS[@]} )); then
    echo "ERROR: task ${tid} produced tmap_idx=${tmap_idx}, but n_maps=${#TMAPS[@]}" >&2
    exit 1
fi

tmap="${MAPS[$tmap_idx]}"
seed="${rep_id}"

mkdir -p "${DEMO_DIR}"
log_msg "GONE2 task
TASK_ID=${tid}
TMAP=${tmap} (tmap_idx=${tmap_idx}/${#TMAPS[@]})
SEED=${seed}
POPS=${POPS[*]}"

for pop in "${POPS[@]}"; do
    vcf="${VCF_DIR}/${tmap}_${seed}_${pop}.vcf.gz"
    vcf_in="${VCF_DIR}/${tmap}_${seed}_${pop}.vcf"
    gone2_out="${DEMO_DIR}/${tmap}_${seed}_${pop}"

    # uncompress vcf file
    bcftools view \
        -i 'F_MISSING <= 0.9' \
        -Ov -o "$vcf_in" \
        "$vcf"

    # run gone2. note, the ADMIX cohort runs with -x flag
    if [[ "${pop}" == *ADMIX* ]]; then
        /home1/karatas/software/GONE2/gone2 \
            -r "$FLAT_RATE" \
            -t "$THREADS" \
            -o "$gone2_out" \
            -S "$seed" \
            -x \
            "$vcf_in"
    else
        /home1/karatas/software/GONE2/gone2 \
            -r "$FLAT_RATE" \
            -t "$THREADS" \
            -o "$gone2_out" \
            -S "$seed" \
            -x \
            "$vcf_in"
    fi

    # delete uncompressed file
    rm -rf "$vcf_in"
done