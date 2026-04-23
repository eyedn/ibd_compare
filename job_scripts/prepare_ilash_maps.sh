#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           prepare_ilash_maps.sh
###############################################################################

#SBATCH --cpus-per-task=1
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
ml conda
source /apps/conda/miniforge3/25.3.0/etc/profile.d/conda.sh
conda activate stdpopsim
export PATH="/home1/karatas/.conda/envs/stdpopsim/bin:$PATH"
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh

: "${SLURM_ARRAY_TASK_ID:?ERROR: run as a Slurm array}"

# inputs
ped_dir="$1"
ilash_map_dir="$2"
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
    reference_map="${ref_dir}/HapMapII_GRCh38/plink.chr${chr}.GRCh38.UCSC.format.map"
elif [[ "$imap" == "DeCodeSexAveraged_GRCh38" ]]; then
    reference_map="${ref_dir}/deCODE/addendum/sex-averaged_gmaps/sex_av_chr${chr}_hg38.gmap"
elif [[ "$imap" == "Flat1Mb1cM" ]]; then
    reference_map="${ref_dir}/GRCh38.p14/recomb_maps_1Mb1cM/chr${chr}.map"
else
    echo "ERROR: unknown imap=${imap}" >&2
    exit 1
fi

ped_map="${ped_dir}/${tmap}_${seed}_${pop}.map"
[[ -s "${ped_map}" ]] || { echo "ERROR: missing PED map ${ped_map}" >&2; exit 1; }
[[ -s "${reference_map}" ]] || {
    echo "ERROR: missing reference map ${reference_map}" >&2
    exit 1
}

map_dir="${ilash_map_dir}/${tmap}/seed${seed}/${pop}/imap${imap}"
mkdir -p "${map_dir}"
output_map="${map_dir}/chr${chr}.ilash.map"

if [[ -s "${output_map}" ]]; then
    log_msg "ilash map exists (${output_map}); skipping"
    exit 0
fi

log_msg "prepare ilash maps
task_id=${tid}
tmap=${tmap}
imap=${imap}
seed=${seed}
pop=${pop}
reference_map=${reference_map}
ped_map=${ped_map}
output_map=${output_map}"

python /home1/karatas/ibd_compare/job_scripts/build_ilash_map.py \
    "${reference_map}" \
    "${ped_map}" \
    "${output_map}"

[[ -s "${output_map}" ]] || {
    echo "ERROR: expected ilash map missing: ${output_map}" >&2
    exit 1
}

log_msg "prepared ilash map ${output_map}"
