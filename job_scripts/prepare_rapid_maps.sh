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

# inputs
vcf_dir="$1"
rapid_map_dir="$2"
ref_dir="$3"
chr="$4"
num_reps="$5"
shift 5

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
    int_genetic_map="${ref_dir}/RaPID_format_maps/HapMapII_GRCh38/plink.chr${chr}.GRCh38.UCSC.reformated.map"
elif [[ "$imap" == "DeCodeSexAveraged_GRCh38" ]]; then
    int_genetic_map="${ref_dir}/RaPID_format_maps/deCODE/chr${chr}.reformated.map"
elif [[ "$imap" == "Flat1Mb1cM" ]]; then
    int_genetic_map="${ref_dir}/RaPID_format_maps/1Mb1cM/chr${chr}.reformated.map"
else
    echo "ERROR: unknown imap=${imap}" >&2
    exit 1
fi

[[ -s "${int_genetic_map}" ]] || {
    echo "ERROR: missing intermediate rapid map ${int_genetic_map}" >&2
    exit 1
}

vcf="${vcf_dir}/${tmap}_${seed}_${pop}.vcf.gz"
[[ -s "${vcf}" ]] || { echo "ERROR: missing VCF ${vcf}" >&2; exit 1; }

map_dir="${rapid_map_dir}/${tmap}/seed${seed}/${pop}/imap${imap}"
mkdir -p "${map_dir}"

filtered_map="${map_dir}/chr${chr}.RaPIDfiltered.map"
final_map="${map_dir}/chr${chr}.RaPIDfinal.map"

if [[ -s "${final_map}" ]]; then
    log_msg "rapid final map exists (${final_map}); skipping"
    exit 0
fi

log_msg "prepare rapid maps
task_id=${tid}
tmap=${tmap}
imap=${imap}
seed=${seed}
pop=${pop}
vcf=${vcf}
int_genetic_map=${int_genetic_map}
filtered_map=${filtered_map}
final_map=${final_map}"

python ~/software/RaPID/filter_mapping_file.py \
    "${int_genetic_map}" \
    "${filtered_map}"
python ~/software/RaPID/interpolate_loci.py \
    "${filtered_map}" \
    "${vcf}" \
    "${final_map}"

[[ -s "${final_map}" ]] || {
    echo "ERROR: expected rapid final map missing: ${final_map}" >&2
    exit 1
}

log_msg "prepared rapid map ${final_map}"
