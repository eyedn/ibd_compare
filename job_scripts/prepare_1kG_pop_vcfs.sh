#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           prepare_1kG_pop_vcfs.sh
###############################################################################

# SBATCH commands if appropriate that were not specified in another script that directly calls this file

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
module purge
ml gcc/13.3.0 htslib/1.19.1 bcftools/1.19
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh

: "${SLURM_ARRAY_TASK_ID:?ERROR: run as a Slurm array}"
threads="${SLURM_CPUS_PER_TASK:-1}"

# parse inputs
source_vcf_dir="$1"
allsamps_dir="$2"
onekg_vcf_dir="$3"
shift 3
chrs=( "$@" )

[[ ${#chrs[@]} -gt 0 ]] || { echo "ERROR: no chromosomes provided" >&2; exit 1; }

tid="${SLURM_ARRAY_TASK_ID}"
chr_idx=$(( tid - 1 ))

if (( chr_idx < 0 || chr_idx >= ${#chrs[@]} )); then
    echo "ERROR: task ${tid} produced chr_idx=${chr_idx}" >&2
    exit 1
fi

chr="${chrs[$chr_idx]}"
chr_dir="${onekg_vcf_dir}/chr${chr}"
mkdir -p "${chr_dir}"

# find exactly one chr-specific source VCF
matches=()
while IFS= read -r source_vcf; do
    base_name="$(basename "${source_vcf}")"
    if [[ "${base_name}" =~ (^|[^[:digit:]])chr${chr}([^[:digit:]]|$) ]]; then
        matches+=( "${source_vcf}" )
    fi
done < <(find "${source_vcf_dir}" -maxdepth 1 -type f -name "*.vcf.gz" | sort)

if (( ${#matches[@]} == 0 )); then
    echo "ERROR: no source VCF matched chr${chr} in ${source_vcf_dir}" >&2
    exit 1
fi
if (( ${#matches[@]} > 1 )); then
    printf 'ERROR: multiple source VCFs matched chr%s:\n' "${chr}" >&2
    printf '  %s\n' "${matches[@]}" >&2
    exit 1
fi

source_vcf="${matches[0]}"
source_tbi="${source_vcf}.tbi"
[[ -s "${source_tbi}" ]] || {
    echo "ERROR: missing source VCF index ${source_tbi}" >&2
    exit 1
}

pops=(ASW CEU YRI)
log_msg "prepare_1kG_pop_vcfs task
task_id=${tid}
chr=${chr}
source_vcf=${source_vcf}
source_tbi=${source_tbi}
allsamps_dir=${allsamps_dir}
out_dir=${chr_dir}
pops=${pops[*]}"

for pop in "${pops[@]}"; do
    sample_file="${allsamps_dir}/${pop}_hg38.txt"
    out_vcf="${chr_dir}/${pop}.vcf.gz"
    out_tbi="${out_vcf}.tbi"

    [[ -s "${sample_file}" ]] || {
        echo "ERROR: missing sample file ${sample_file}" >&2
        exit 1
    }

    if [[ -s "${out_vcf}" && -s "${out_tbi}" ]]; then
        log_msg "subset VCF reused (${out_vcf})"
        continue
    fi

    log_msg "subsetting chr${chr} ${pop}
sample_file=${sample_file}
out_vcf=${out_vcf}"
    bcftools view \
        -S "${sample_file}" \
        --threads "${threads}" \
        -Oz \
        -o "${out_vcf}" \
        "${source_vcf}"
    tabix -f -p vcf "${out_vcf}"

    [[ -s "${out_vcf}" && -s "${out_tbi}" ]] || {
        echo "ERROR: expected subset outputs missing for ${out_vcf}" >&2
        exit 1
    }
done

log_msg "prepare_1kG_pop_vcfs done for chr${chr}"
