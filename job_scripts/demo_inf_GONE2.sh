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
threads="${SLURM_CPUS_PER_TASK:-1}"

# inputs
vcf_dir="$1"
demo_dir="$2"
chr="$3"
flat_rate="$4"
num_reps="$5"
shift 5

# read tmaps until "--", then read pops after
tmaps=()
while [[ $# -gt 0 && "$1" != "--" ]]; do
    tmaps+=( "$1" )
    shift
done
if [[ $# -eq 0 ]]; then
    echo "ERROR: expected '--' separator between tmaps and pops" >&2
    exit 1
fi
shift 1
pops=( "$@" )

if [[ ${#tmaps[@]} -eq 0 ]]; then
    echo "ERROR: no tmaps provided" >&2
    exit 1
fi
if [[ ${#pops[@]} -eq 0 ]]; then
    echo "ERROR: no pops provided" >&2
    exit 1
fi

tid="${SLURM_ARRAY_TASK_ID}"
tmap_idx=$(( (tid - 1) / num_reps ))
rep_id=$(( (tid - 1) % num_reps + 1 ))

if (( tmap_idx < 0 || tmap_idx >= ${#tmaps[@]} )); then
    echo "ERROR: task ${tid} produced tmap_idx=${tmap_idx}, but n_maps=${#tmaps[@]}" >&2
    exit 1
fi

tmap="${tmaps[$tmap_idx]}"
seed="${rep_id}"

common_gone2_dir="${demo_dir}/gone2"
mkdir -p "${common_gone2_dir}"
log_msg "GONE2 task
task_id=${tid}
tmap=${tmap} (tmap_idx=${tmap_idx}/${#tmaps[@]})
seed=${seed}
chr=${chr}
flat_rate=${flat_rate}
pops=${pops[*]}"

for pop in "${pops[@]}"; do
    log_msg "detected pop=${pop}"
    vcf="${vcf_dir}/${tmap}_${seed}_${pop}.vcf.gz"
    vcf_in="${vcf_dir}/${tmap}_${seed}_${pop}_tmp.vcf"
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
            -r "$flat_rate" \
            -t "$threads" \
            -o "${gone2_out}_x" \
            -S "$seed" \
            -x \
            "$vcf_in"
    fi

    log_msg "running gone without -x"
    /home1/karatas/software/GONE2/gone2 \
        -r "$flat_rate" \
        -t "$threads" \
        -o "$gone2_out" \
        -S "$seed" \
        "$vcf_in"

    # delete uncompressed file
    rm -rf "$vcf_in"
done