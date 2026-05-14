#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           prepare_1kG_rapid_maps.sh
###############################################################################

# SBATCH commands if appropriate that were not specified in another script that directly calls this file

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

# parse inputs
vcf_dir="$1"
rapid_map_dir="$2"
ref_dir="$3"
shift 3

chrs=()
while [[ $# -gt 0 && "$1" != "--" ]]; do
    chrs+=( "$1" )
    shift
done
[[ $# -gt 0 && "$1" == "--" ]] || {
    echo "ERROR: expected '--' after chromosomes" >&2
    exit 1
}
shift

imaps=()
while [[ $# -gt 0 && "$1" != "--" ]]; do
    imaps+=( "$1" )
    shift
done
[[ $# -gt 0 && "$1" == "--" ]] || {
    echo "ERROR: expected '--' after imaps" >&2
    exit 1
}
shift

pops=( "$@" )

[[ ${#chrs[@]} -gt 0 ]] || { echo "ERROR: no chromosomes provided" >&2; exit 1; }
[[ ${#imaps[@]} -gt 0 ]] || { echo "ERROR: no imaps provided" >&2; exit 1; }
[[ ${#pops[@]} -gt 0 ]] || { echo "ERROR: no pops provided" >&2; exit 1; }

# decode task id -> (chr, imap, pop)
tid="${SLURM_ARRAY_TASK_ID}"
tid0=$(( tid - 1 ))

pop_idx=$(( tid0 % ${#pops[@]} ))
tmp_id=$(( tid0 / ${#pops[@]} ))
imap_idx=$(( tmp_id % ${#imaps[@]} ))
chr_idx=$(( tmp_id / ${#imaps[@]} ))

if (( chr_idx < 0 || chr_idx >= ${#chrs[@]} )); then
    echo "ERROR: task ${tid} produced chr_idx=${chr_idx}" >&2
    exit 1
fi
if (( imap_idx < 0 || imap_idx >= ${#imaps[@]} )); then
    echo "ERROR: task ${tid} produced imap_idx=${imap_idx}" >&2
    exit 1
fi
if (( pop_idx < 0 || pop_idx >= ${#pops[@]} )); then
    echo "ERROR: task ${tid} produced pop_idx=${pop_idx}" >&2
    exit 1
fi

chr="${chrs[$chr_idx]}"
imap="${imaps[$imap_idx]}"
pop="${pops[$pop_idx]}"

# choose the same rapid reference-map source used by simulation prep
if [[ "${imap}" == "HapMapII_GRCh38" ]]; then
    rapid_int_map="${ref_dir}/RaPID_format_maps/HapMapII_GRCh38/plink.chr${chr}.GRCh38.UCSC.reformated.map"
elif [[ "${imap}" == "DeCodeSexAveraged_GRCh38" ]]; then
    rapid_int_map="${ref_dir}/RaPID_format_maps/deCODE/chr${chr}.reformated.map"
elif [[ "${imap}" == "Flat1Mb1cM" ]]; then
    rapid_int_map="${ref_dir}/RaPID_format_maps/1Mb1cM/chr${chr}.reformated.map"
else
    echo "ERROR: unknown imap=${imap}" >&2
    exit 1
fi

[[ -s "${rapid_int_map}" ]] || {
    echo "ERROR: missing intermediate rapid map ${rapid_int_map}" >&2
    exit 1
}

vcf="${vcf_dir}/chr${chr}/${pop}.vcf.gz"
[[ -s "${vcf}" ]] || { echo "ERROR: missing subset VCF ${vcf}" >&2; exit 1; }

mkdir -p "${rapid_map_dir}"
rapid_final_map="${rapid_map_dir}/pop${pop}_chr${chr}_imap${imap}.map"

if [[ -s "${rapid_final_map}" ]]; then
    log_msg "rapid map reused (${rapid_final_map})"
    exit 0
fi

tmp_root="${TMPDIR:-/tmp}"
tmp_dir="$(mktemp -d \
    "${tmp_root}/prepare_1kG_rapid_maps.${pop}_chr${chr}_${imap}.XXXXXX")"
cleanup() {
    rm -rf "${tmp_dir}"
}
trap cleanup EXIT

rapid_filtered_map="${tmp_dir}/pop${pop}_chr${chr}_imap${imap}.filteredIntermediate"

log_msg "prepare_1kG_rapid_maps task
task_id=${tid}
chr=${chr}
imap=${imap}
pop=${pop}
vcf=${vcf}
rapid_int_map=${rapid_int_map}
rapid_final_map=${rapid_final_map}"

python ~/software/RaPID/filter_mapping_file.py \
    "${rapid_int_map}" \
    "${rapid_filtered_map}"
python ~/software/RaPID/interpolate_loci.py \
    "${rapid_filtered_map}" \
    "${vcf}" \
    "${rapid_final_map}"

[[ -s "${rapid_final_map}" ]] || {
    echo "ERROR: expected rapid final map missing: ${rapid_final_map}" >&2
    exit 1
}

log_msg "prepared rapid map ${rapid_final_map}"
