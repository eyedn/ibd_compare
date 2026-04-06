#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           filter_format_hapibd.sh
###############################################################################

#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=8G
#SBATCH --time=1-00:00:00
#SBATCH --partition=qcb
#SBATCH --account=jazlynmo_738
#SBATCH --nodes=1
#SBATCH --output=/home1/karatas/logs/hapFF.%A_%a.out
#SBATCH --error=/home1/karatas/logs/hapFF.%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=karatas@usc.edu

set -euo pipefail
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh

: "${SLURM_ARRAY_TASK_ID:?ERROR: run as a Slurm array}"

# inputs
OUTDIR="$1"
VCF_DIR="$2"
COMMON_FF_DIR="$3"
REF_DIR="$4"
CHR="$5"
MIN_BP_LEN="$6"
MIN_CM_LEN="$7"
HOT_BIN_SIZE="$8"
HOT_BIN_SD="$9"
NUM_REPS="${10}"

# hapibd params (from submit script)
HAP_MIN_SEED="${11}"
HAP_MAX_GAP="${12}"
HAP_MIN_EXTEND="${13}"
HAP_MIN_OUTPUT="${14}"
HAP_MIN_MARKERS="${15}"
HAP_MIN_MAC="${16}"
shift 16

# read TMAPs until "--", then IMAPS until "--", then POPS
TMAPS=()
while [[ $# -gt 0 && "$1" != "--" ]]; do TMAPS+=( "$1" ); shift; done
[[ $# -gt 0 && "$1" == "--" ]] || { echo "ERROR: expected '--' after TMAPS" >&2; exit 1; }
shift

IMAPS=()
while [[ $# -gt 0 && "$1" != "--" ]]; do IMAPS+=( "$1" ); shift; done
[[ $# -gt 0 && "$1" == "--" ]] || { echo "ERROR: expected '--' after IMAPS" >&2; exit 1; }
shift

POPS=( "$@" )

[[ ${#TMAPS[@]} -gt 0 ]] || { echo "ERROR: no TMAPS provided" >&2; exit 1; }
[[ ${#IMAPS[@]} -gt 0 ]] || { echo "ERROR: no IMAPS provided" >&2; exit 1; }
[[ ${#POPS[@]} -gt 0 ]] || { echo "ERROR: no POPS provided" >&2; exit 1; }

# decode array task -> (tmap, imap, seed)
tid="${SLURM_ARRAY_TASK_ID}"
tid0=$(( tid - 1 ))

seed_id=$(( tid0 % NUM_REPS + 1 ))
tmp_id=$(( tid0 / NUM_REPS ))

imap_idx=$(( tmp_id % ${#IMAPS[@]} ))
tmap_idx=$(( tmp_id / ${#IMAPS[@]} ))

if (( tmap_idx < 0 || tmap_idx >= ${#TMAPS[@]} )); then
    echo "ERROR: task ${tid} produced tmap_idx=${tmap_idx}, n_tmaps=${#TMAPS[@]}" >&2
    exit 1
fi
if (( imap_idx < 0 || imap_idx >= ${#IMAPS[@]} )); then
    echo "ERROR: task ${tid} produced imap_idx=${imap_idx}, n_imaps=${#IMAPS[@]}" >&2
    exit 1
fi

tmap="${TMAPS[$tmap_idx]}"
imap="${IMAPS[$imap_idx]}"
seed="${seed_id}"

mkdir -p "$COMMON_FF_DIR"

# standardized labels
if [[ "$tmap" == *HapMapII_GRCh38* ]]; then tmap_label="hapmap"
elif [[ "$tmap" == *DeCodeSexAveraged_GRCh38* ]]; then tmap_label="decode"
elif [[ "$tmap" == *DefaultUniform* ]]; then tmap_label="flat"
else tmap_label="$tmap"
fi

if [[ "$imap" == *HapMapII_GRCh38* ]]; then
    imap_label="hapmap"
    hapibd_gmap="${REF_DIR}/HapMapII_GRCh38/plink.chr${CHR}.GRCh38.UCSC.format.map"
elif [[ "$imap" == *DeCodeSexAveraged_GRCh38* ]]; then
    imap_label="decode"
    hapibd_gmap="${REF_DIR}/deCODE/addendum/sex-averaged_gmaps/sex_av_chr${CHR}_hg38.gmap"
elif [[ "$imap" == *Flat1Mb1cM* ]]; then
    imap_label="1Mb1cM"
    hapibd_gmap="${REF_DIR}/GRCh38.p14/recomb_maps_1Mb1cM/chr${CHR}.map"
else
    imap_label="$imap"
    hapibd_gmap=""
fi

[[ -n "$hapibd_gmap" && -s "$hapibd_gmap" ]] || { echo "ERROR: missing hapibd_gmap: $hapibd_gmap" >&2; exit 1; }

# hapibd dir base (use const params)
hapibd_dir_base="hapibd_sim_tmap${tmap}_imap${imap}_seed${HAP_MIN_SEED}_gap${HAP_MAX_GAP}_extend${HAP_MIN_EXTEND}_out${HAP_MIN_OUTPUT}_mark${HAP_MIN_MARKERS}_mac${HAP_MIN_MAC}"
hapibd_dir="${OUTDIR}/${hapibd_dir_base}"

for pop in "${POPS[@]}"; do
    vcf="${VCF_DIR}/${tmap}_${seed}_${pop}.vcf.gz"

    hapibd_nonbed="${hapibd_dir}/${pop}_chr${CHR}_seed${seed}.ibd.gz"
    hapibd_ibdends_prefix="${hapibd_dir}/${pop}_chr${CHR}_seed${seed}.ibdends"
    hapibd_ibdends_ibd="${hapibd_ibdends_prefix}_ibd.ibdends.bed.gz"

    # (your existing bed made elsewhere)
    hapibd_ibd="${hapibd_dir}/${pop}_chr${CHR}_seed${seed}_ibd.bed.gz"

    hapibd_ff_dir="${COMMON_FF_DIR}/${hapibd_dir_base}"
    hapibd_out="${hapibd_ff_dir}/hapibd_${tmap_label}_${imap_label}_${pop}_${seed}.bed"
    hapibd_hot_out="${hapibd_ff_dir}/hapibd_${tmap_label}_${imap_label}_${pop}_${seed}.hot.bed"
    hapibd_ibdends_out="${hapibd_ff_dir}/hapibd_${tmap_label}_${imap_label}_${pop}_${seed}.ibdends.bed"
    hapibd_ibdendshot_out="${hapibd_ff_dir}/hapibd_${tmap_label}_${imap_label}_${pop}_${seed}.ibdendshot.bed"

    [[ -s "$vcf" ]] || { echo "ERROR: missing vcf: $vcf" >&2; exit 1; }
    [[ -s "$hapibd_nonbed" ]] || { echo "ERROR: missing hapibd ibd.gz: $hapibd_nonbed" >&2; exit 1; }
    [[ -s "$hapibd_ibd" ]] || { echo "ERROR: missing hapibd bed.gz: $hapibd_ibd" >&2; exit 1; }

    mkdir -p "$hapibd_ff_dir"

    # hapibd -> ibd-ends
    module purge
    ml legacy/CentOS7 jdk/17.0.5
    log_msg "hapibd: filtering segments using ibd-ends (pop=${pop} seed=${seed} tmap=${tmap_label} imap=${imap_label})"
    java -Xmx46G -jar ~/software/ibd-ends/ibd-ends.jar \
        gt="$vcf" \
        ibd="$hapibd_nonbed" \
        map="$hapibd_gmap" \
        nthreads="${SLURM_CPUS_PER_TASK:-1}" \
        out="$hapibd_ibdends_prefix"

    log_msg "hapibd: generating ${hapibd_ibdends_ibd} from ${hapibd_ibdends_prefix}.ibd.gz"
    zcat "${hapibd_ibdends_prefix}.ibd.gz" | \
        awk 'BEGIN{OFS="\t"} NR>1 {print $5,$10,$11,$1,$2,$3,$4,$12}' | \
        gzip -c > "$hapibd_ibdends_ibd"

    # format/filter (R)
    module purge
    ml rstats/4.5.1

    log_msg "hapibd: format+filter raw calls"
    Rscript "/home1/karatas/ibd_compare/job_scripts/filter_and_format.R" \
        "$MIN_BP_LEN" "$MIN_CM_LEN" "$HOT_BIN_SIZE" "$HOT_BIN_SD" \
        "$hapibd_ibd" "$hapibd_gmap" \
        "$hapibd_out" "$hapibd_hot_out" \
        "TRUE"
    gzip -f "$hapibd_out" "$hapibd_hot_out"

    log_msg "hapibd: format+filter ibd-ends adjusted calls"
    Rscript "/home1/karatas/ibd_compare/job_scripts/filter_and_format.R" \
        "$MIN_BP_LEN" "$MIN_CM_LEN" "$HOT_BIN_SIZE" "$HOT_BIN_SD" \
        "$hapibd_ibdends_ibd" "$hapibd_gmap" \
        "$hapibd_ibdends_out" "$hapibd_ibdendshot_out" \
        "TRUE"
    gzip -f "$hapibd_ibdends_out" "$hapibd_ibdendshot_out"
done

log_msg "hapibd pipeline done"
