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
#SBATCH --output=/home1/karatas/logs/tvHap.%A_%a.%x.out
#SBATCH --error=/home1/karatas/logs/tvHap.%A_%a.%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=karatas@usc.edu

set -euo pipefail
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh
: "${SLURM_ARRAY_TASK_ID:?ERROR: run as a Slurm array}"

# inputs
COMMON_FF_DIR="$1"
POP_INFO_DIR="$2"
TVI_DIR="$3"
CHR="$4"
declare -i SAMPLE_SIZE="$5"
NUM_REPS="$6"

# hapibd params (from submit script)
HAP_MIN_SEED="$7"
HAP_MAX_GAP="$8"
HAP_MIN_EXTEND="$9"
HAP_MIN_OUTPUT="${10}"
HAP_MIN_MARKERS="${11}"
HAP_MIN_MAC="${12}"
shift 12

# read TMAPS until "--", then IMAPS until "--", then FILTERING until "--", then POPS
TMAPS=()
while [[ $# -gt 0 && "$1" != "--" ]]; do TMAPS+=( "$1" ); shift; done
[[ $# -gt 0 && "$1" == "--" ]] || { echo "ERROR: expected '--' after TMAPS" >&2; exit 1; }
shift

IMAPS=()
while [[ $# -gt 0 && "$1" != "--" ]]; do IMAPS+=( "$1" ); shift; done
[[ $# -gt 0 && "$1" == "--" ]] || { echo "ERROR: expected '--' after IMAPS" >&2; exit 1; }
shift

FILTERING=()
while [[ $# -gt 0 && "$1" != "--" ]]; do FILTERING+=( "$1" ); shift; done
[[ $# -gt 0 && "$1" == "--" ]] || { echo "ERROR: expected '--' after FILTERING" >&2; exit 1; }
shift

POPS=( "$@" )

# decode task id -> pop + filt + imap + tmap (NO software dimension)
tid="${SLURM_ARRAY_TASK_ID}"
tid0=$(( tid - 1 ))

pop_idx=$(( tid0 % ${#POPS[@]} ))
tmp_id=$(( tid0 / ${#POPS[@]} ))

filt_idx=$(( tmp_id % ${#FILTERING[@]} ))
tmp_id=$(( tmp_id / ${#FILTERING[@]} ))

imap_idx=$(( tmp_id % ${#IMAPS[@]} ))
tmap_idx=$(( tmp_id / ${#IMAPS[@]} ))

# bounds checks
if (( pop_idx < 0 || pop_idx >= ${#POPS[@]} )); then echo "ERROR pop_idx=$pop_idx" >&2; exit 1; fi
if (( filt_idx < 0 || filt_idx >= ${#FILTERING[@]} )); then echo "ERROR filt_idx=$filt_idx" >&2; exit 1; fi
if (( imap_idx < 0 || imap_idx >= ${#IMAPS[@]} )); then echo "ERROR imap_idx=$imap_idx" >&2; exit 1; fi
if (( tmap_idx < 0 || tmap_idx >= ${#TMAPS[@]} )); then echo "ERROR tmap_idx=$tmap_idx" >&2; exit 1; fi

pop="${POPS[$pop_idx]}"
tmap="${TMAPS[$tmap_idx]}"
imap="${IMAPS[$imap_idx]}"
filt="${FILTERING[$filt_idx]}"
software="hapibd"

# labels (unchanged)
if [[ "$tmap" == *HapMapII_GRCh38* ]]; then tmap_label="hapmap"
elif [[ "$tmap" == *DeCodeSexAveraged_GRCh38* ]]; then tmap_label="decode"
elif [[ "$tmap" == *DefaultUniform* ]]; then tmap_label="flat"
else tmap_label="$tmap"
fi

if [[ "$imap" == *HapMapII_GRCh38* ]]; then imap_label="hapmap"
elif [[ "$imap" == *DeCodeSexAveraged_GRCh38* ]]; then imap_label="decode"
elif [[ "$imap" == *Flat1Mb1cM* ]]; then imap_label="1Mb1cM"
else imap_label="$imap"
fi

software_dir_base="hapibd_sim_tmap${tmap}_imap${imap}_seed${HAP_MIN_SEED}_gap${HAP_MAX_GAP}_extend${HAP_MIN_EXTEND}_out${HAP_MIN_OUTPUT}_mark${HAP_MIN_MARKERS}_mac${HAP_MIN_MAC}"
mkdir -p "${TVI_DIR}/${software_dir_base}"

log_msg "detected...
tmap=${tmap}
imap=${imap}
software=${software}
filt=${filt}
pop=${pop}
NUM_REPS=${NUM_REPS}
software_dir_base=${software_dir_base}"

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
    "$manifest" "$tmap_label" "$CHR" "$SAMPLE_SIZE" "$agg_out" "$filt" "$pop"
