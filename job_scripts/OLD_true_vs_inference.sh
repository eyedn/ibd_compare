#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           true_vs_inference.sh
###############################################################################

#SBATCH -p qcb
#SBATCH -A jazlynmo_738
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=1-00:00:00
#SBATCH --output=/home1/karatas/logs/tviIBD.%A_%a.%x.out
#SBATCH --error=/home1/karatas/logs/tviIBD.%A_%a.%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=karatas@usc.edu

set -euo pipefail
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh

# constants
: "${SLURM_ARRAY_TASK_ID:?ERROR: run as a Slurm array}"

# hapibd params (to find correct ibd file)
MIN_SEED=1.0
MAX_GAP=2000
MIN_EXTEND=1.0
MIN_OUTPUT=2.0
MIN_MARKERS=100
MIN_MAC=2

# rapid params (to find correct ibd file)
MIN_LEN=2
WINDOW_SIZE=100
NUM_RUNS=10
NUM_SUCCESSES=2

# inputs
COMMON_FF_DIR="$1"
POP_INFO_DIR="$2"
TVI_DIR="$3"
CHR="$4"
declare -i SAMPLE_SIZE="$5"
NUM_REPS="$6"
shift 6

# read TMAPS until "--", then IMAPS until "--", then SOFTWARE until "--",
# then FILTERING until "--", then POPS (rest)
TMAPS=()
while [[ $# -gt 0 && "$1" != "--" ]]; do
    TMAPS+=( "$1" )
    shift
done
[[ $# -gt 0 && "$1" == "--" ]] || { echo "ERROR: expected '--' after TMAPS" >&2; exit 1; }
shift

IMAPS=()
while [[ $# -gt 0 && "$1" != "--" ]]; do
    IMAPS+=( "$1" )
    shift
done
[[ $# -gt 0 && "$1" == "--" ]] || { echo "ERROR: expected '--' after IMAPS" >&2; exit 1; }
shift

SOFTWARE=()
while [[ $# -gt 0 && "$1" != "--" ]]; do
    SOFTWARE+=( "$1" )
    shift
done
[[ $# -gt 0 && "$1" == "--" ]] || { echo "ERROR: expected '--' after SOFTWARE" >&2; exit 1; }
shift

FILTERING=()
while [[ $# -gt 0 && "$1" != "--" ]]; do
    FILTERING+=( "$1" )
    shift
done
[[ $# -gt 0 && "$1" == "--" ]] || { echo "ERROR: expected '--' after FILTERING" >&2; exit 1; }
shift

POPS=( "$@" )

# decode task id to which pop + tmap + imap + software + filtering
tid="${SLURM_ARRAY_TASK_ID}"
tid0=$(( tid - 1 ))

pop_idx=$(( tid0 % ${#POPS[@]} ))
tmp_id=$(( tid0 / ${#POPS[@]} ))

filt_idx=$(( tmp_id % ${#FILTERING[@]} ))
tmp_id=$(( tmp_id / ${#FILTERING[@]} ))

software_idx=$(( tmp_id % ${#SOFTWARE[@]} ))
tmp_id=$(( tmp_id / ${#SOFTWARE[@]} ))

imap_idx=$(( tmp_id % ${#IMAPS[@]} ))
tmap_idx=$(( tmp_id / ${#IMAPS[@]} ))

if (( pop_idx < 0 || pop_idx >= ${#POPS[@]} )); then
    echo "ERROR: task ${tid} produced pop_idx=${pop_idx}, but n_pops=${#POPS[@]}" >&2
    exit 1
fi

if (( tmap_idx < 0 || tmap_idx >= ${#TMAPS[@]} )); then
    echo "ERROR: task ${tid} produced tmap_idx=${tmap_idx}, but n_tmaps=${#TMAPS[@]}" >&2
    exit 1
fi
if (( imap_idx < 0 || imap_idx >= ${#IMAPS[@]} )); then
    echo "ERROR: task ${tid} produced imap_idx=${imap_idx}, but n_imaps=${#IMAPS[@]}" >&2
    exit 1
fi
if (( software_idx < 0 || software_idx >= ${#SOFTWARE[@]} )); then
    echo "ERROR: task ${tid} produced software_idx=${software_idx}, but n_software=${#SOFTWARE[@]}" >&2
    exit 1
fi
if (( filt_idx < 0 || filt_idx >= ${#FILTERING[@]} )); then
    echo "ERROR: task ${tid} produced filt_idx=${filt_idx}, but n_filtering=${#FILTERING[@]}" >&2
    exit 1
fi

pop="${POPS[$pop_idx]}"
tmap="${TMAPS[$tmap_idx]}"
imap="${IMAPS[$imap_idx]}"
software="${SOFTWARE[$software_idx]}"
filt="${FILTERING[$filt_idx]}"

log_msg "detected...
tmap=${tmap}
imap=${imap}
software=${software}
filt=${filt}
pop=${pop}
NUM_REPS=${NUM_REPS}"

# convert to standardized tmap label
if [[ "$tmap" == *HapMapII_GRCh38* ]]; then 
    tmap_label="hapmap"
elif [[ "$tmap" == *DeCodeSexAveraged_GRCh38* ]]; then
    tmap_label="decode"
elif [[ "$tmap" == *DefaultUniform* ]]; then
    tmap_label="flat"
fi

# convert to standardized imap label
if [[ "$imap" == *HapMapII_GRCh38* ]]; then 
    imap_label="hapmap"
elif [[ "$imap" == *DeCodeSexAveraged_GRCh38* ]]; then
    imap_label="decode"
elif [[ "$imap" == *Flat1Mb1cM* ]]; then
    imap_label="1Mb1cM"
fi

# build software_dir_base
if [[ "$software" == "hapibd" ]]; then
    software_dir_base="hapibd_sim_tmap${tmap}_imap${imap}_seed${MIN_SEED}_gap${MAX_GAP}_extend${MIN_EXTEND}_out${MIN_OUTPUT}_mark${MIN_MARKERS}_mac${MIN_MAC}"
else
    software_dir_base="rapid_sim_tmap${tmap}_imap${imap}_d${MIN_LEN}_w${WINDOW_SIZE}_r${NUM_RUNS}_s${NUM_SUCCESSES}"
fi

# loop through POP and SEED to process each inference vs. true ibd segments
mkdir -p "${TVI_DIR}/${software_dir_base}"

# manifest listing one row per seed
manifest="${TVI_DIR}/${software_dir_base}/${software}_${tmap_label}_${imap_label}_${filt}_${pop}.manifest.tsv"
printf "seed\ttrue_ibd\tinf_ibd\tcoverage_out\tunrels\n" > "$manifest"

for seed in $(seq 1 "$NUM_REPS"); do
    true_ibd="${COMMON_FF_DIR}/true_ibd/true_${tmap_label}_${tmap_label}_${pop}_${seed}.bed.gz"
    unrels="${POP_INFO_DIR}/${tmap}_${seed}_${pop}_unrel_samples.txt"

    if [[ "$filt" == "none" ]]; then
    inf_ibd="${COMMON_FF_DIR}/${software_dir_base}/${software}_${tmap_label}_${imap_label}_${pop}_${seed}.bed.gz"
    else
    inf_ibd="${COMMON_FF_DIR}/${software_dir_base}/${software}_${tmap_label}_${imap_label}_${pop}_${seed}.${filt}.bed.gz"
    fi

    coverage_out="${TVI_DIR}/${software_dir_base}/${software}_${tmap_label}_${imap_label}_${filt}_${pop}_${seed}.bed.gz"

    [[ -s "$true_ibd" ]] || { log_msg "ERROR: missing true_ibd=$true_ibd"; exit 1; }
    [[ -s "$inf_ibd"  ]] || { log_msg "ERROR: missing inf_ibd=$inf_ibd"; exit 1; }
    [[ -s "$unrels"   ]] || { log_msg "ERROR: missing unrels=$unrels"; exit 1; }

    log_msg "seed=${seed} | bedtools coverage"

    module purge
    ml gcc/13.3.0 bedtools2/2.31.1

    bedtools coverage -sorted \
    -a <(
        zcat "$true_ibd" | tail -n +2 | \
        awk 'BEGIN{OFS="\t"} $8 != $10 { key=$1 "|" $8 "_" $9 "|" $10 "_" $11; $1=key; print }' | \
        LC_ALL=C sort -k1,1 -k2,2n -k3,3n
    ) \
    -b <(
        zcat "$inf_ibd" | tail -n +2 | \
        awk 'BEGIN{OFS="\t"} $8 != $10 { key=$1 "|" $8 "_" $9 "|" $10 "_" $11; $1=key; print }' | \
        LC_ALL=C sort -k1,1 -k2,2n -k3,3n
    ) | \
    cut -f1-13,15 | \
    awk 'BEGIN{OFS="\t"} { split($1, a, /\|/); $1 = a[1]; print }' | \
    gzip -c > "$coverage_out"

    log_msg "wrote $coverage_out"

    # record this replicate in the manifest
    printf "%s\t%s\t%s\t%s\t%s\n" \
    "$seed" "$true_ibd" "$inf_ibd" "$coverage_out" "$unrels" >> "$manifest"
done

log_msg "finished bedtools for pop=${pop}; manifest=$manifest"
agg_out="${TVI_DIR}/${software_dir_base}/${software}_${tmap_label}_${imap_label}_${filt}_${pop}.agg.tsv"

log_msg "running true_vs_inference_aggregate.R"
module purge
ml rstats/4.5.1
Rscript "/home1/karatas/ibd_compare/job_scripts/true_vs_inference_aggregate.R" \
    "$manifest" \
    "$tmap_label" \
    "$CHR" \
    "$SAMPLE_SIZE" \
    "$agg_out" \
    "$filt" \
    "$pop"
