#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           inference_1kG_hapibd.sh
###############################################################################

# SBATCH commands if appropriate that were not specified in another script that directly calls this file

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --time=1-00:00:00
#SBATCH --partition=qcb
#SBATCH --account=jazlynmo_738
#SBATCH --nodes=1
#SBATCH --output=/home1/karatas/logs/hapibd.%A_%a.out
#SBATCH --error=/home1/karatas/logs/hapibd.%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=karatas@usc.edu

set -euo pipefail
module purge
ml openjdk/25.0.2
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh

: "${SLURM_ARRAY_TASK_ID:?ERROR: run as a Slurm array}"

# parse inputs
min_seed="$1"
max_gap="$2"
min_extend="$3"
min_output="$4"
min_markers="$5"
min_mac="$6"
outdir="$7"
vcf_dir="$8"
ref_dir="$9"
shift 9

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

# choose the same hap-IBD map source used by simulation inference
if [[ "${imap}" == "HapMapII_GRCh38" ]]; then
    genetic_map="${ref_dir}/HapMapII_GRCh38/plink.chr${chr}.GRCh38.UCSC.format.map"
elif [[ "${imap}" == "DeCodeSexAveraged_GRCh38" ]]; then
    genetic_map="${ref_dir}/deCODE/addendum/sex-averaged_gmaps/sex_av_chr${chr}_hg38.gmap"
elif [[ "${imap}" == "Flat1Mb1cM" ]]; then
    genetic_map="${ref_dir}/GRCh38.p14/recomb_maps_1Mb1cM/chr${chr}.map"
else
    echo "ERROR: unknown imap=${imap}" >&2
    exit 1
fi

[[ -s "${genetic_map}" ]] || {
    echo "ERROR: missing hap-IBD map ${genetic_map}" >&2
    exit 1
}

ibd_dir="${outdir}/hapibd_1kG_imap${imap}_seed${min_seed}_gap${max_gap}_extend${min_extend}_out${min_output}_mark${min_markers}_mac${min_mac}"
mkdir -p "${ibd_dir}"

log_msg "inference_1kG_hapibd task
task_id=${tid}
chr=${chr}
imap=${imap}
genetic_map=${genetic_map}
ibd_dir=${ibd_dir}
pops=${pops[*]}"

for pop in "${pops[@]}"; do
    vcf="${vcf_dir}/chr${chr}/${pop}.vcf.gz"
    prefix="${ibd_dir}/${pop}_chr${chr}"
    ibd_gz="${prefix}.ibd.gz"
    ibd_bed="${prefix}_ibd.bed.gz"

    [[ -s "${vcf}" ]] || { echo "ERROR: missing subset VCF ${vcf}" >&2; exit 1; }

    if [[ -s "${ibd_gz}" && -s "${ibd_bed}" ]]; then
        log_msg "hap-IBD outputs reused (${ibd_bed})"
        continue
    fi

    if [[ ! -s "${ibd_gz}" ]]; then
        log_msg "calling hap-IBD
pop=${pop}
chr=${chr}
imap=${imap}
gt=${vcf}
map=${genetic_map}
out=${prefix}"

        rm -f \
            "${prefix}.log" \
            "${prefix}.hbd.gz" \
            "${prefix}.ibd.gz" \
            "${prefix}.ibd.gz.tbi" \
            2>/dev/null || true
        java -Xmx14G -jar ~/software/hap-ibd/hap-ibd.jar \
            min-seed="${min_seed}" \
            max-gap="${max_gap}" \
            min-extend="${min_extend}" \
            min-output="${min_output}" \
            min-markers="${min_markers}" \
            min-mac="${min_mac}" \
            nthreads="${SLURM_CPUS_PER_TASK:-1}" \
            gt="${vcf}" \
            map="${genetic_map}" \
            out="${prefix}"
    fi

    [[ -s "${ibd_gz}" ]] || {
        echo "ERROR: expected hap-IBD output missing: ${ibd_gz}" >&2
        exit 1
    }

    tmpdir="${SCRATCH:-/tmp}/bedtmp_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID:-$$}"
    mkdir -p "${tmpdir}"
    tmp_bed="${ibd_bed}.tmp"

    if [[ $(zcat "${ibd_gz}" | wc -l) -eq 0 ]]; then
        gzip -c /dev/null > "${tmp_bed}"
    else
        zcat "${ibd_gz}" | \
        awk -F'\t' 'BEGIN{OFS="\t"} {print $5,$6,$7,$1,$2,$3,$4,$8}' | \
        sort -k1,1 -k2,2n -k3,3n -T "${tmpdir}" -S 2G | \
        gzip > "${tmp_bed}"
    fi

    mv -f "${tmp_bed}" "${ibd_bed}"
    rm -rf "${tmpdir}"
done

log_msg "called 1kG hap-IBD for chr=${chr}, imap=${imap}"
