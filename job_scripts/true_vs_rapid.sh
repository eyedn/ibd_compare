#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           true_vs_rapid.sh
###############################################################################

#SBATCH -p qcb
#SBATCH -A jazlynmo_738
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=1-00:00:00
#SBATCH --output=/home1/karatas/logs/tvRap.%A_%a.%x.out
#SBATCH --error=/home1/karatas/logs/tvRap.%A_%a.%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=karatas@usc.edu

set -euo pipefail
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh
: "${SLURM_ARRAY_TASK_ID:?ERROR: run as a Slurm array}"

# inputs
common_ff_dir="$1"
pop_info_dir="$2"
tvi_dir="$3"
chr="$4"
declare -i sample_size="$5"
num_reps="$6"

# rapid params (from submit script)
rap_min_len="$7"
rap_window_size="$8"
rap_num_runs="$9"
rap_num_successes="${10}"
shift 10

# read tmaps until "--", then imaps until "--", then filtering until "--", then pops
tmaps=()
while [[ $# -gt 0 && "$1" != "--" ]]; do tmaps+=( "$1" ); shift; done
[[ $# -gt 0 && "$1" == "--" ]] || { echo "ERROR: expected '--' after tmaps" >&2; exit 1; }
shift

imaps=()
while [[ $# -gt 0 && "$1" != "--" ]]; do imaps+=( "$1" ); shift; done
[[ $# -gt 0 && "$1" == "--" ]] || { echo "ERROR: expected '--' after imaps" >&2; exit 1; }
shift

filtering=()
while [[ $# -gt 0 && "$1" != "--" ]]; do filtering+=( "$1" ); shift; done
[[ $# -gt 0 && "$1" == "--" ]] || { echo "ERROR: expected '--' after filtering" >&2; exit 1; }
shift

pops=( "$@" )

# decode task id -> pop + filt + imap + tmap
tid="${SLURM_ARRAY_TASK_ID}"
tid0=$(( tid - 1 ))

pop_idx=$(( tid0 % ${#pops[@]} ))
tmp_id=$(( tid0 / ${#pops[@]} ))

filt_idx=$(( tmp_id % ${#filtering[@]} ))
tmp_id=$(( tmp_id / ${#filtering[@]} ))

imap_idx=$(( tmp_id % ${#imaps[@]} ))
tmap_idx=$(( tmp_id / ${#imaps[@]} ))

if (( pop_idx < 0 || pop_idx >= ${#pops[@]} )); then echo "ERROR pop_idx=$pop_idx" >&2; exit 1; fi
if (( filt_idx < 0 || filt_idx >= ${#filtering[@]} )); then echo "ERROR filt_idx=$filt_idx" >&2; exit 1; fi
if (( imap_idx < 0 || imap_idx >= ${#imaps[@]} )); then echo "ERROR imap_idx=$imap_idx" >&2; exit 1; fi
if (( tmap_idx < 0 || tmap_idx >= ${#tmaps[@]} )); then echo "ERROR tmap_idx=$tmap_idx" >&2; exit 1; fi

pop="${pops[$pop_idx]}"
tmap="${tmaps[$tmap_idx]}"
imap="${imaps[$imap_idx]}"
filt="${filtering[$filt_idx]}"
software="rapid"

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

software_dir_base="rapid_sim_tmap${tmap}_imap${imap}_d${rap_min_len}_w${rap_window_size}_r${rap_num_runs}_s${rap_num_successes}"
mkdir -p "${tvi_dir}/${software_dir_base}"

log_msg "detected...
tmap=${tmap}
imap=${imap}
software=${software}
filt=${filt}
pop=${pop}
num_reps=${num_reps}
software_dir_base=${software_dir_base}"

manifest="${tvi_dir}/${software_dir_base}/${software}_${tmap_label}_${imap_label}_${filt}_${pop}.manifest.tsv"
printf "seed\ttrue_ibd\tinf_ibd\tcoverage_out\tunrels\n" > "$manifest"

for seed in $(seq 1 "$num_reps"); do
    true_ibd="${common_ff_dir}/true_ibd/true_${tmap_label}_${tmap_label}_${pop}_${seed}.bed.gz"
    unrels="${pop_info_dir}/${tmap}_${seed}_${pop}_unrel_samples.txt"

    if [[ "$filt" == "none" ]]; then
        inf_ibd="${common_ff_dir}/${software_dir_base}/${software}_${tmap_label}_${imap_label}_${pop}_${seed}.bed.gz"
    else
        inf_ibd="${common_ff_dir}/${software_dir_base}/${software}_${tmap_label}_${imap_label}_${pop}_${seed}.${filt}.bed.gz"
    fi

    coverage_out="${tvi_dir}/${software_dir_base}/${software}_${tmap_label}_${imap_label}_${filt}_${pop}_${seed}.bed.gz"

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
        lc_all=C sort -k1,1 -k2,2n -k3,3n
        ) \
        -b <(
        zcat "$inf_ibd" | tail -n +2 | \
        awk 'BEGIN{OFS="\t"} $8 != $10 { key=$1 "|" $8 "_" $9 "|" $10 "_" $11; $1=key; print }' | \
        lc_all=C sort -k1,1 -k2,2n -k3,3n
        ) | \
        cut -f1-13,15 | \
        awk 'BEGIN{OFS="\t"} { split($1, a, /\|/); $1 = a[1]; print }' | \
        gzip -c > "$coverage_out"

    # record this replicate in the manifest
    printf "%s\t%s\t%s\t%s\t%s\n" \
        "$seed" "$true_ibd" "$inf_ibd" "$coverage_out" "$unrels" >> "$manifest"
done

log_msg "finished bedtools for pop=${pop}; manifest=$manifest"
agg_out="${tvi_dir}/${software_dir_base}/${software}_${tmap_label}_${imap_label}_${filt}_${pop}.agg.tsv"

log_msg "running true_vs_inference_aggregate.R"
module purge
ml rstats/4.5.1
Rscript "/home1/karatas/ibd_compare/job_scripts/true_vs_inference_aggregate.R" \
    "$manifest" "$tmap_label" "$chr" "$sample_size" "$agg_out" "$filt" "$pop"
