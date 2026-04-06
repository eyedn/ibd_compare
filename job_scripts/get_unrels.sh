#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           get_unrels.sh
###############################################################################

#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=1-00:00:00
#SBATCH --partition=qcb
#SBATCH --account=jazlynmo_738
#SBATCH --nodes=1
#SBATCH --output=/home1/karatas/logs/slurm.%A_%a.%x.out
#SBATCH --error=/home1/karatas/logs/slurm.%A_%a.%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=karatas@usc.edu

set -euo pipefail
ml gcc/13.3.0 bcftools/1.19 htslib/1.19.1 plink2/2.00a4.3
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh

: "${SLURM_ARRAY_TASK_ID:?ERROR: run as a Slurm array}"
THREADS="${SLURM_CPUS_PER_TASK:-1}"

# inputs
KIN_CUTOFF="$1"
VCF_DIR="$2"
POP_INFO_DIR="$3"
NUM_REPS="$4"
shift 4

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

# for each pop, use its pop-specific vcf to find unrelated individuals
plink_out_dir="${POP_INFO_DIR}/plink2_king_out"
mkdir -p "${plink_out_dir}"
log_msg "get_unrels task
TASK_ID=${tid}
GMAP=${gmap}
SEED=${seed}
KIN_CUTOFF=${KIN_CUTOFF}
POPS=${POPS[*]}"
for pop in "${POPS[@]}"; do
  in_vcf="${VCF_DIR}/${gmap}_${seed}_${pop}.vcf.gz"
  [[ -s "${in_vcf}" ]] || { echo "ERROR: missing input VCF ${in_vcf}" >&2; exit 1; }

  out_prefix="${plink_out_dir}/${gmap}_${seed}_${pop}"
  king_in="${out_prefix}.king.cutoff.in.id"

  unrel_ids="${POP_INFO_DIR}/${gmap}_${seed}_${pop}_unrel_samples.txt"
  unrel_vcf="${VCF_DIR}/${gmap}_${seed}_${pop}_unrels.vcf.gz"

  if [[ -s "${unrel_ids}" && -s "${unrel_vcf}" ]]; then
    log_msg "unrels exist (${unrel_vcf}); skipping ${gmap}_${seed}_${pop}"
    continue
  fi

  log_msg "running plink2 KING cutoff=${KIN_CUTOFF} on ${gmap}_${seed}_${pop}"
  plink2 \
    --vcf "${in_vcf}" \
    --king-cutoff "${KIN_CUTOFF}" \
    --make-king \
    --make-king-table \
    --threads "${THREADS}" \
    --out "${out_prefix}"

  [[ -s "${king_in}" ]] || { echo "ERROR: expected ${king_in} not found" >&2; exit 1; }

  # extract IIDs; handle header ("#IID") safely
awk 'NF>0 && $1!="IID" && $1!="#IID" {print $1}' "${king_in}" | sort -u > "${unrel_ids}"
  log_msg "wrote unrel IDs -> ${unrel_ids} ($(wc -l < "${unrel_ids}") samples)"

  log_msg "creating unrels VCF for ${gmap}_${seed}_${pop}"
  bcftools view \
    -S "${unrel_ids}" \
    --threads "${THREADS}" \
    -Oz -o "${unrel_vcf}" \
    "${in_vcf}"
  tabix -p vcf "${unrel_vcf}"
  log_msg "created ${unrel_vcf}"
done

log_msg "created unrels for all pops for ${gmap} seed=${seed}"
