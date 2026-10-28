#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           inference_1kG_rapid.sh
###############################################################################


#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --time=1-00:00:00
#SBATCH --partition=qcb
#SBATCH --account=jazlynmo_738
#SBATCH --nodes=1
#SBATCH --output=/home1/karatas/logs/1kGrapid/%A_%a.%x.out
#SBATCH --error=/home1/karatas/logs/1kGrapid/%A_%a.%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=karatas@usc.edu

set -euo pipefail
ml purge
ml gcc/13.3.0 htslib/1.19.1
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh

: "${SLURM_ARRAY_TASK_ID:?ERROR: run as a Slurm array}"

# parse inputs
min_len="$1"
window_size="$2"
num_runs="$3"
num_successes="$4"
vcf_dir="$5"
rapid_map_dir="$6"
outdir="$7"
shift 7

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

# decode task id -> (chr, imap)
tid="${SLURM_ARRAY_TASK_ID}"
tid0=$(( tid - 1 ))
imap_idx=$(( tid0 % ${#imaps[@]} ))
chr_idx=$(( tid0 / ${#imaps[@]} ))

if (( chr_idx < 0 || chr_idx >= ${#chrs[@]} )); then
    echo "ERROR: task ${tid} produced chr_idx=${chr_idx}" >&2
    exit 1
fi
if (( imap_idx < 0 || imap_idx >= ${#imaps[@]} )); then
    echo "ERROR: task ${tid} produced imap_idx=${imap_idx}" >&2
    exit 1
fi

chr="${chrs[$chr_idx]}"
imap="${imaps[$imap_idx]}"

ibd_dir="${outdir}/rapid_1kG_imap${imap}_d${min_len}_w${window_size}_r${num_runs}_s${num_successes}"
mkdir -p "${ibd_dir}"

log_msg "inference_1kG_rapid task
task_id=${tid}
chr=${chr}
imap=${imap}
ibd_dir=${ibd_dir}
pops=${pops[*]}"

for pop in "${pops[@]}"; do
    vcf="${vcf_dir}/chr${chr}/${pop}.vcf.gz"
    genetic_map="${rapid_map_dir}/pop${pop}_chr${chr}_imap${imap}.map"
    ibd_out="${ibd_dir}/${pop}_chr${chr}_out"
    rapid_result="${ibd_out}/results.max.gz"
    ibd_bed="${ibd_dir}/${pop}_chr${chr}_ibd.bed.gz"

    [[ -s "${vcf}" ]] || { echo "ERROR: missing subset VCF ${vcf}" >&2; exit 1; }
    [[ -s "${genetic_map}" ]] || {
        echo "ERROR: missing rapid map ${genetic_map}" >&2
        exit 1
    }

    if [[ -s "${rapid_result}" && -s "${ibd_bed}" ]]; then
        log_msg "rapid outputs reused (${ibd_bed})"
        continue
    fi

    if [[ ! -s "${rapid_result}" ]]; then
        log_msg "calling rapid
pop=${pop}
chr=${chr}
imap=${imap}
input_file_vcf_compressed=${vcf}
genetic_mapping_file=${genetic_map}
output_folder=${ibd_out}"

        rm -rf "${ibd_out}"
        "${HOME}/software/RaPID/RaPID_v.1.7" \
            -i "${vcf}" \
            -g "${genetic_map}" \
            -d "${min_len}" \
            -w "${window_size}" \
            -r "${num_runs}" \
            -s "${num_successes}" \
            -o "${ibd_out}"
    fi

    [[ -s "${rapid_result}" ]] || {
        echo "ERROR: expected rapid output missing: ${rapid_result}" >&2
        exit 1
    }

    tmpdir="/home1/karatas/scratch/bedtmp_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID:-$$}"
    mkdir -p "${tmpdir}"
    tmp_bed="${ibd_bed}.tmp"

    if [[ $(zcat "${rapid_result}" | wc -l) -eq 0 ]]; then
        gzip -c /dev/null > "${tmp_bed}"
    else
        zcat "${rapid_result}" | \
        awk -F'\t' 'BEGIN{OFS="\t"}
            $6 ~ /^-?[0-9]+$/ && $7 ~ /^-?[0-9]+$/ && $6 >= 0 && $7 >= $6 {
            print $1,$6,$7,$2,$4,$3,$5,$8,$9,$10
            }' | \
        sort -k1,1 -k2,2n -k3,3n -T "${tmpdir}" -S 2G | \
        gzip > "${tmp_bed}"
    fi

    mv -f "${tmp_bed}" "${ibd_bed}"
    rm -rf "${tmpdir}"
done

log_msg "called 1kG rapid for chr=${chr}, imap=${imap}"
