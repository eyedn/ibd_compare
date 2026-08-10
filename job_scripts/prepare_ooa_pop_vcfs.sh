#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           prepare_ooa_pop_vcfs.sh
###############################################################################

#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=1-00:00:00
#SBATCH --partition=qcb
#SBATCH --account=jazlynmo_738
#SBATCH --nodes=1
#SBATCH --output=/home1/karatas/logs/prepIn/%A_%a.%x.out
#SBATCH --error=/home1/karatas/logs/prepIn/%A_%a.%x.err

set -euo pipefail
module purge
ml gcc/13.3.0 htslib/1.19.1 bcftools/1.19
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh

: "${SLURM_ARRAY_TASK_ID:?ERROR: run as a Slurm array}"
threads="${SLURM_CPUS_PER_TASK:-1}"
ooa_vcf_dir="$1"
vcf_dir="$2"
pop_info_dir="$3"
chr="$4"
sample_size="$5"
num_reps="$6"
contract_file="$7"
id_contract="$8"
shift 8

tmaps=()
while [[ $# -gt 0 && "$1" != "--" ]]; do
    tmaps+=("$1")
    shift
done
[[ $# -gt 0 && "$1" == "--" ]] || {
    echo "ERROR: expected '--' separator between maps and populations" >&2
    exit 1
}
shift
pops=("$@")
[[ ${#tmaps[@]} -gt 0 ]] || { echo "ERROR: no maps provided" >&2; exit 1; }
[[ ${#pops[@]} -gt 0 ]] || { echo "ERROR: no populations provided" >&2; exit 1; }

tid="${SLURM_ARRAY_TASK_ID}"
tmap_idx=$(( (tid - 1) / num_reps ))
rep=$(( (tid - 1) % num_reps + 1 ))
if (( tmap_idx < 0 || tmap_idx >= ${#tmaps[@]} )); then
    echo "ERROR: task ${tid} produced tmap_idx=${tmap_idx}" >&2
    exit 1
fi
tmap="${tmaps[$tmap_idx]}"
source_vcf="${ooa_vcf_dir}/${tmap}_${rep}_chr${chr}_all.vcf.gz"
[[ -s "${source_vcf}" ]] || { echo "ERROR: missing OOA VCF ${source_vcf}" >&2; exit 1; }
[[ -s "${source_vcf}.tbi" || -s "${source_vcf}.csi" ]] || {
    echo "ERROR: missing OOA VCF index for ${source_vcf}" >&2
    exit 1
}

mkdir -p "${vcf_dir}" "${pop_info_dir}"
log_msg "prepare OOA population VCFs: tmap=${tmap}; rep=${rep}; chr=${chr}; source_vcf=${source_vcf}; pops=${pops[*]}"
for pop in "${pops[@]}"; do
    sample_ids="${pop_info_dir}/${tmap}_${rep}_chr${chr}_${pop}_samples.txt"
    out_vcf="${vcf_dir}/${tmap}_${rep}_chr${chr}_${pop}.vcf.gz"
    if [[ ! -s "${sample_ids}" ]]; then
        : > "${sample_ids}"
        for sample_num in $(seq 1 "${sample_size}"); do
            printf '%s_%s\n' "${pop}" "${sample_num}" >> "${sample_ids}"
        done
    fi
    expected_ids="${sample_ids}.expected"
    {
        for sample_num in $(seq 1 "${sample_size}"); do
            printf '%s_%s\n' "${pop}" "${sample_num}"
        done
    } > "${expected_ids}"
    cmp -s "${sample_ids}" "${expected_ids}" || {
        echo "ERROR: stale or mixed sample IDs in ${sample_ids}" >&2
        exit 1
    }
    rm -f "${expected_ids}"
    if [[ ! -s "${out_vcf}" || ( ! -s "${out_vcf}.tbi" && ! -s "${out_vcf}.csi" ) ]]; then
        bcftools view -S "${sample_ids}" --threads "${threads}" \
            -Oz -o "${out_vcf}" "${source_vcf}"
        tabix -p vcf "${out_vcf}"
    fi
    bcftools query -l "${out_vcf}" > "${out_vcf}.samples"
    cmp -s "${sample_ids}" "${out_vcf}.samples" || {
        echo "ERROR: VCF ID contract mismatch in ${out_vcf}" >&2
        exit 1
    }
    rm -f "${out_vcf}.samples"
done
printf '%s\n' "${id_contract}" > "${contract_file}.tmp"
mv -f "${contract_file}.tmp" "${contract_file}"
