#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           prepare_rapid_ilash_maps.sh
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

# inputs
vcf_dir="$1"
ped_dir="$2"
rapid_map_dir="$3"
ilash_map_dir="$4"
ref_dir="$5"
chr="$6"
num_reps="$7"
shift 7

tmaps=()
while [[ $# -gt 0 && "$1" != "--" ]]; do
    tmaps+=( "$1" )
    shift
done
[[ $# -gt 0 && "$1" == "--" ]] || {
    echo "ERROR: expected '--' after tmaps" >&2
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

[[ ${#tmaps[@]} -gt 0 ]] || { echo "ERROR: no tmaps provided" >&2; exit 1; }
[[ ${#imaps[@]} -gt 0 ]] || { echo "ERROR: no imaps provided" >&2; exit 1; }
[[ ${#pops[@]} -gt 0 ]] || { echo "ERROR: no pops provided" >&2; exit 1; }

# decode array task -> (tmap, imap, seed, pop)
tid="${SLURM_ARRAY_TASK_ID}"
tid0=$(( tid - 1 ))

pop_idx=$(( tid0 % ${#pops[@]} ))
tmp_id=$(( tid0 / ${#pops[@]} ))

seed_id=$(( tmp_id % num_reps + 1 ))
tmp_id=$(( tmp_id / num_reps ))

imap_idx=$(( tmp_id % ${#imaps[@]} ))
tmap_idx=$(( tmp_id / ${#imaps[@]} ))

if (( pop_idx < 0 || pop_idx >= ${#pops[@]} )); then
    echo "ERROR: task ${tid} produced pop_idx=${pop_idx}" >&2
    exit 1
fi
if (( imap_idx < 0 || imap_idx >= ${#imaps[@]} )); then
    echo "ERROR: task ${tid} produced imap_idx=${imap_idx}" >&2
    exit 1
fi
if (( tmap_idx < 0 || tmap_idx >= ${#tmaps[@]} )); then
    echo "ERROR: task ${tid} produced tmap_idx=${tmap_idx}" >&2
    exit 1
fi

tmap="${tmaps[$tmap_idx]}"
imap="${imaps[$imap_idx]}"
seed="${seed_id}"
pop="${pops[$pop_idx]}"

if [[ "$imap" == "HapMapII_GRCh38" ]]; then
    rapid_int_map="${ref_dir}/RaPID_format_maps/HapMapII_GRCh38/plink.chr${chr}.GRCh38.UCSC.reformated.map"
    ilash_reference_map="${ref_dir}/HapMapII_GRCh38/plink.chr${chr}.GRCh38.UCSC.format.map"
elif [[ "$imap" == "DeCodeSexAveraged_GRCh38" ]]; then
    rapid_int_map="${ref_dir}/RaPID_format_maps/deCODE/chr${chr}.reformated.map"
    ilash_reference_map="${ref_dir}/deCODE/addendum/sex-averaged_gmaps/sex_av_chr${chr}_hg38.gmap"
elif [[ "$imap" == "Flat1Mb1cM" ]]; then
    rapid_int_map="${ref_dir}/RaPID_format_maps/1Mb1cM/chr${chr}.reformated.map"
    ilash_reference_map="${ref_dir}/GRCh38.p14/recomb_maps_1Mb1cM/chr${chr}.map"
else
    echo "ERROR: unknown imap=${imap}" >&2
    exit 1
fi

[[ -s "${rapid_int_map}" ]] || {
    echo "ERROR: missing intermediate rapid map ${rapid_int_map}" >&2
    exit 1
}
[[ -s "${ilash_reference_map}" ]] || {
    echo "ERROR: missing ilash reference map ${ilash_reference_map}" >&2
    exit 1
}

vcf="${vcf_dir}/${tmap}_${seed}_${pop}.vcf.gz"
ped_map="${ped_dir}/${tmap}_${seed}_${pop}.map"
[[ -s "${vcf}" ]] || { echo "ERROR: missing VCF ${vcf}" >&2; exit 1; }
[[ -s "${ped_map}" ]] || { echo "ERROR: missing PED map ${ped_map}" >&2; exit 1; }

rapid_dir="${rapid_map_dir}/${tmap}/seed${seed}/${pop}/imap${imap}"
ilash_dir="${ilash_map_dir}/${tmap}/seed${seed}/${pop}/imap${imap}"
mkdir -p "${rapid_dir}" "${ilash_dir}"

rapid_filtered_map="${rapid_dir}/chr${chr}.RaPIDfiltered.map"
rapid_final_map="${rapid_dir}/chr${chr}.RaPIDfinal.map"
ilash_output_map="${ilash_dir}/chr${chr}.ilash.map"

need_rapid=true
need_ilash=true

if [[ -s "${rapid_final_map}" ]]; then
    need_rapid=false
fi
if [[ -s "${ilash_output_map}" ]]; then
    need_ilash=false
fi

if [[ "${need_rapid}" == "false" && "${need_ilash}" == "false" ]]; then
    log_msg "rapid and ilash maps exist for ${tmap} seed=${seed} ${pop} ${imap}; skipping"
    exit 0
fi

log_msg "prepare rapid+ilash maps
task_id=${tid}
tmap=${tmap}
imap=${imap}
seed=${seed}
pop=${pop}
vcf=${vcf}
ped_map=${ped_map}
rapid_int_map=${rapid_int_map}
rapid_final_map=${rapid_final_map}
ilash_reference_map=${ilash_reference_map}
ilash_output_map=${ilash_output_map}"

if [[ "${need_rapid}" == "true" ]]; then
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
else
    log_msg "rapid map reused (${rapid_final_map})"
fi

if [[ "${need_ilash}" == "true" ]]; then
    python /home1/karatas/ibd_compare/job_scripts/build_ilash_map.py \
        "${ilash_reference_map}" \
        "${ped_map}" \
        "${ilash_output_map}"
    [[ -s "${ilash_output_map}" ]] || {
        echo "ERROR: expected ilash map missing: ${ilash_output_map}" >&2
        exit 1
    }
    log_msg "prepared ilash map ${ilash_output_map}"
else
    log_msg "ilash map reused (${ilash_output_map})"
fi
