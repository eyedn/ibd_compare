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
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --partition=qcb
#SBATCH --account=jazlynmo_738
#SBATCH --nodes=1
#SBATCH --output=/home1/karatas/logs/gone2.%A_%a.%x.out
#SBATCH --error=/home1/karatas/logs/gone2.%A_%a.%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=karatas@usc.edu

set -euo pipefail
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh

ml purge
ml gcc/13.3.0 htslib/1.19.1 bcftools/1.19

: "${SLURM_ARRAY_TASK_ID:?ERROR: run as a Slurm array}"
THREADS="${SLURM_CPUS_PER_TASK:-1}"

# inputs
VCF_DIR="$1"
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

tmap="${TMAPS[$tmap_idx]}"
seed="${rep_id}"

common_gone2_dir="${DEMO_DIR}/gone2"
mkdir -p "${common_gone2_dir}"
log_msg "GONE2 task
TASK_ID=${tid}
TMAP=${tmap} (tmap_idx=${tmap_idx}/${#TMAPS[@]})
SEED=${seed}
CHR=${CHR}
FLAT_RATE=${FLAT_RATE}
POPS=${POPS[*]}"

for pop in "${POPS[@]}"; do
    log_msg "detected pop=${pop}"
    vcf="${VCF_DIR}/${tmap}_${seed}_${pop}.vcf.gz"
    vcf_in="${VCF_DIR}/${tmap}_${seed}_${pop}_tmp.vcf"
    gone2_out="${common_gone2_dir}/${tmap}_${seed}_${pop}"

    # uncompress vcf file
    bcftools view \
        -i 'F_MISSING <= 0.9' \
        -Ov -o "$vcf_in" \
        "$vcf"

    snp_count=$(grep -v "^#" "$vcf_in" | wc -l)
    log_msg "runing gone2 on ${vcf_in}
    snp_count=${snp_count}"

    # run gone2. note, the ADMIX cohort runs with -x flag
    if [[ "$pop" == *ADMIX* ]]; then
        log_msg "running gone with -x"
        /home1/karatas/software/GONE2/gone2 \
            -r "$FLAT_RATE" \
            -t "$THREADS" \
            -o "${gone2_out}_x" \
            -S "$seed" \
            -x \
            "$vcf_in"
    fi

    log_msg "running gone without -x"
    /home1/karatas/software/GONE2/gone2 \
        -r "$FLAT_RATE" \
        -t "$THREADS" \
        -o "$gone2_out" \
        -S "$seed" \
        "$vcf_in"

    # delete uncompressed file
    rm -rf "$vcf_in"
done