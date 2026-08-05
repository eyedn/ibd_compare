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
threads="${SLURM_CPUS_PER_TASK:-1}"

# inputs
kin_cutoff="$1"
vcf_dir="$2"
pop_info_dir="$3"
chr="$4"
num_reps="$5"
shift 5

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

all_outputs_exist=true
for pop in "${pops[@]}"; do
  unrel_ids="${pop_info_dir}/${gmap}_${rep}_chr${chr}_${pop}_unrel_samples.txt"
  unrel_vcf="${vcf_dir}/${gmap}_${rep}_chr${chr}_${pop}_unrels.vcf.gz"
  if [[ ! -s "${unrel_ids}" || ! -s "${unrel_vcf}" ]]; then
    all_outputs_exist=false
    break
  fi
done

if [[ "${all_outputs_exist}" == "true" ]]; then
  log_msg "all unrel outputs exist for ${gmap} rep=${rep}; skipping"
  exit 0
fi

# for each pop, use its pop-specific vcf to find unrelated individuals
plink_out_dir="${pop_info_dir}/plink2_king_out"
mkdir -p "${plink_out_dir}"
log_msg "get_unrels task
task_id=${tid}
gmap=${gmap}
rep=${rep}
kin_cutoff=${kin_cutoff}
pops=${pops[*]}"
for pop in "${pops[@]}"; do
  in_vcf="${vcf_dir}/${gmap}_${rep}_chr${chr}_${pop}.vcf.gz"
  [[ -s "${in_vcf}" ]] || { echo "ERROR: missing input VCF ${in_vcf}" >&2; exit 1; }

  out_prefix="${plink_out_dir}/${gmap}_${rep}_chr${chr}_${pop}"
  king_in="${out_prefix}.king.cutoff.in.id"

  unrel_ids="${pop_info_dir}/${gmap}_${rep}_chr${chr}_${pop}_unrel_samples.txt"
  unrel_vcf="${vcf_dir}/${gmap}_${rep}_chr${chr}_${pop}_unrels.vcf.gz"

  if [[ -s "${unrel_ids}" && -s "${unrel_vcf}" ]]; then
    log_msg "unrels exist (${unrel_vcf}); skipping ${gmap}_${rep}_${pop}"
    continue
  fi

  log_msg "running plink2 KING cutoff=${kin_cutoff} on ${gmap}_${rep}_${pop}"
  plink2 \
    --vcf "${in_vcf}" \
    --king-cutoff "${kin_cutoff}" \
    --make-king \
    --make-king-table \
    --threads "${threads}" \
    --out "${out_prefix}"

  [[ -s "${king_in}" ]] || { echo "ERROR: expected ${king_in} not found" >&2; exit 1; }

  # extract IIDs; handle header ("#IID") safely
awk 'NF>0 && $1!="IID" && $1!="#IID" {print $1}' "${king_in}" | sort -u > "${unrel_ids}"
  log_msg "wrote unrel IDs -> ${unrel_ids} ($(wc -l < "${unrel_ids}") samples)"

  log_msg "creating unrels VCF for ${gmap}_${rep}_${pop}"
  bcftools view \
    -S "${unrel_ids}" \
    --threads "${threads}" \
    -Oz -o "${unrel_vcf}" \
    "${in_vcf}"
  tabix -p vcf "${unrel_vcf}"
  log_msg "created ${unrel_vcf}"
done

log_msg "created unrels for all pops for ${gmap} rep=${rep}"
