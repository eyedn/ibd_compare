#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           prepare_rapid_maps.sh
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
ml gcc/13.3.0 htslib/1.19.1 bcftools/1.19 conda
source /apps/conda/miniforge3/25.3.0/etc/profile.d/conda.sh
conda activate rpgg_proj
export PATH="/home1/karatas/.conda/envs/rpgg_proj/bin:$PATH"
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh

: "${SLURM_ARRAY_TASK_ID:?ERROR: run as a Slurm array}"
vcf_dir="$1"
rapid_map_dir="$2"
ref_dir="$3"
chr="$4"
num_reps="$5"
contract_file="$6"
id_contract="$7"
shift 7
source /home1/karatas/ibd_compare/other_scripts/require_sim_id_contract.sh \
    "$contract_file" "$id_contract"

tmaps=()
while [[ $# -gt 0 && "$1" != "--" ]]; do
    tmaps+=("$1")
    shift
done
[[ $# -gt 0 && "$1" == "--" ]] || { echo "ERROR: expected '--' after maps" >&2; exit 1; }
shift
imaps=()
while [[ $# -gt 0 && "$1" != "--" ]]; do
    imaps+=("$1")
    shift
done
[[ $# -gt 0 && "$1" == "--" ]] || { echo "ERROR: expected '--' after maps" >&2; exit 1; }
shift
pops=("$@")
[[ ${#tmaps[@]} -gt 0 && ${#imaps[@]} -gt 0 && ${#pops[@]} -gt 0 ]] || {
    echo "ERROR: maps and populations are required" >&2
    exit 1
}

tid0=$(( SLURM_ARRAY_TASK_ID - 1 ))
pop_idx=$(( tid0 % ${#pops[@]} ))
tmp_id=$(( tid0 / ${#pops[@]} ))
rep=$(( tmp_id % num_reps + 1 ))
tmp_id=$(( tmp_id / num_reps ))
imap_idx=$(( tmp_id % ${#imaps[@]} ))
tmap_idx=$(( tmp_id / ${#imaps[@]} ))
if (( pop_idx >= ${#pops[@]} || imap_idx >= ${#imaps[@]} || tmap_idx >= ${#tmaps[@]} )); then
    echo "ERROR: task ${SLURM_ARRAY_TASK_ID} is outside the configured array" >&2
    exit 1
fi
tmap="${tmaps[$tmap_idx]}"
imap="${imaps[$imap_idx]}"
pop="${pops[$pop_idx]}"

case "${imap}" in
    HapMapII_GRCh38)
        rapid_int_map="${ref_dir}/RaPID_format_maps/HapMapII_GRCh38/plink.chr${chr}.GRCh38.UCSC.reformated.map"
        ;;
    DeCodeSexAveraged_GRCh38)
        rapid_int_map="${ref_dir}/RaPID_format_maps/deCODE/chr${chr}.reformated.map"
        ;;
    Flat1Mb1cM)
        rapid_int_map="${ref_dir}/RaPID_format_maps/1Mb1cM/chr${chr}.reformated.map"
        ;;
    *) echo "ERROR: unknown imap=${imap}" >&2; exit 1 ;;
esac

vcf="${vcf_dir}/${tmap}_${rep}_chr${chr}_${pop}.vcf.gz"
map_prefix="pop${pop}_chr${chr}_rep${rep}_tmap${tmap}_imap${imap}"
rapid_filtered_map="${rapid_map_dir}/${map_prefix}.filteredIntermediate"
rapid_final_map="${rapid_map_dir}/${map_prefix}.map"
[[ -s "${rapid_int_map}" ]] || { echo "ERROR: missing map ${rapid_int_map}" >&2; exit 1; }
[[ -s "${vcf}" ]] || { echo "ERROR: missing VCF ${vcf}" >&2; exit 1; }
mkdir -p "${rapid_map_dir}"
if [[ ! -s "${rapid_final_map}" ]]; then
    python ~/software/RaPID/filter_mapping_file.py \
        "${rapid_int_map}" "${rapid_filtered_map}"
    python ~/software/RaPID/interpolate_loci.py \
        "${rapid_filtered_map}" "${vcf}" "${rapid_final_map}"
fi
[[ -s "${rapid_final_map}" ]] || { echo "ERROR: missing RaPID map ${rapid_final_map}" >&2; exit 1; }
log_msg "prepared RaPID map: tmap=${tmap}; imap=${imap}; rep=${rep}; pop=${pop}; map=${rapid_final_map}"
