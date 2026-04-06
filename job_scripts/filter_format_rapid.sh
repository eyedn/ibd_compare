#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           filter_format_rapid.sh
###############################################################################

#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=8G
#SBATCH --time=1-00:00:00
#SBATCH --partition=qcb
#SBATCH --account=jazlynmo_738
#SBATCH --nodes=1
#SBATCH --output=/home1/karatas/logs/rapFF.%A_%a.out
#SBATCH --error=/home1/karatas/logs/rapFF.%A_%a.err
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

# rapid params (from submit script)
RAP_MIN_LEN="${11}"
RAP_WINDOW_SIZE="${12}"
RAP_NUM_RUNS="${13}"
RAP_NUM_SUCCESSES="${14}"
shift 14

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

# NOTE: ibd-ends needs a NON-RaPID-style map; you were using hapibd_gmap for this
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

# rapid dir base (use const params)
rapid_dir_base="rapid_sim_tmap${tmap}_imap${imap}_d${RAP_MIN_LEN}_w${RAP_WINDOW_SIZE}_r${RAP_NUM_RUNS}_s${RAP_NUM_SUCCESSES}"
rapid_dir="${OUTDIR}/${rapid_dir_base}"

for pop in "${POPS[@]}"; do
    vcf="${VCF_DIR}/${tmap}_${seed}_${pop}.vcf.gz"

    rapid_nonbed_dir="${rapid_dir}/${pop}_chr${CHR}_seed${seed}_out"
    rapid_nonbed="${rapid_nonbed_dir}/results.max.gz"
    rapid_nonbed_hapibdformat="${rapid_nonbed_dir}/results.max.hapibdformat.gz"
    rapid_ibdends_prefix="${rapid_dir}/${pop}_chr${CHR}_seed${seed}.ibdends"
    rapid_ibdends_ibd="${rapid_ibdends_prefix}_ibd.ibdends.bed.gz"

    # (your existing bed made elsewhere)
    rapid_ibd="${rapid_dir}/${pop}_chr${CHR}_seed${seed}_ibd.bed.gz"

    # rapid map for cM interpolation during formatting
    rapid_gmap_dir="${rapid_dir}/maps"
    rapid_gmap="${rapid_gmap_dir}/chr${CHR}.${pop}.${seed}.RaPIDfinal.map"

    rapid_ff_dir="${COMMON_FF_DIR}/${rapid_dir_base}"
    rapid_out="${rapid_ff_dir}/rapid_${tmap_label}_${imap_label}_${pop}_${seed}.bed"
    rapid_hot_out="${rapid_ff_dir}/rapid_${tmap_label}_${imap_label}_${pop}_${seed}.hot.bed"
    rapid_ibdends_out="${rapid_ff_dir}/rapid_${tmap_label}_${imap_label}_${pop}_${seed}.ibdends.bed"
    rapid_ibdendshot_out="${rapid_ff_dir}/rapid_${tmap_label}_${imap_label}_${pop}_${seed}.ibdendshot.bed"

    [[ -s "$vcf" ]] || { echo "ERROR: missing vcf: $vcf" >&2; exit 1; }
    [[ -s "$rapid_nonbed" ]] || { echo "ERROR: missing rapid results: $rapid_nonbed" >&2; exit 1; }
    [[ -s "$rapid_ibd" ]] || { echo "ERROR: missing rapid bed.gz: $rapid_ibd" >&2; exit 1; }
    [[ -s "$rapid_gmap" ]] || { echo "ERROR: missing rapid gmap: $rapid_gmap" >&2; exit 1; }

    mkdir -p "$rapid_ff_dir"

    # rapid -> ibd-ends (convert to hapibd-like format first)
    module purge
    ml legacy/CentOS7 jdk/17.0.5

    log_msg "rapid: generating ${rapid_nonbed_hapibdformat} from ${rapid_nonbed}"
    zcat "$rapid_nonbed" | \
        awk 'BEGIN{OFS="\t"} ($6 >= 0 && $7 >= 0 && $6 <= $7) { print $2, $4+1, $3, $5+1, $1, $6, $7, $8 }' | \
        gzip -c > "$rapid_nonbed_hapibdformat"

    log_msg "rapid: filtering segments using ibd-ends (pop=${pop} seed=${seed} tmap=${tmap_label} imap=${imap_label})"
    java -Xmx46G -jar ~/software/ibd-ends/ibd-ends.jar \
        gt="$vcf" \
        ibd="$rapid_nonbed_hapibdformat" \
        map="$hapibd_gmap" \
        nthreads="${SLURM_CPUS_PER_TASK:-1}" \
        out="$rapid_ibdends_prefix"

    log_msg "rapid: generating ${rapid_ibdends_ibd} from ${rapid_ibdends_prefix}.ibd.gz"
    zcat "${rapid_ibdends_prefix}.ibd.gz" | \
        awk 'BEGIN{OFS="\t"} NR>1 {print $5,$10,$11,$1,$2,$3,$4,$12}' | \
        gzip -c > "$rapid_ibdends_ibd"

    # format/filter (R)
    module purge
    ml rstats/4.5.1

    log_msg "rapid: format+filter raw calls"
    Rscript "/home1/karatas/ibd_compare/job_scripts/filter_and_format.R" \
        "$MIN_BP_LEN" "$MIN_CM_LEN" "$HOT_BIN_SIZE" "$HOT_BIN_SD" \
        "$rapid_ibd" "$rapid_gmap" \
        "$rapid_out" "$rapid_hot_out" \
        "FALSE"
    gzip -f "$rapid_out" "$rapid_hot_out"

    log_msg "rapid: format+filter ibd-ends adjusted calls (hapibd-style map here)"
    Rscript "/home1/karatas/ibd_compare/job_scripts/filter_and_format.R" \
        "$MIN_BP_LEN" "$MIN_CM_LEN" "$HOT_BIN_SIZE" "$HOT_BIN_SD" \
        "$rapid_ibdends_ibd" "$hapibd_gmap" \
        "$rapid_ibdends_out" "$rapid_ibdendshot_out" \
        "TRUE"
    gzip -f "$rapid_ibdends_out" "$rapid_ibdendshot_out"
done

log_msg "rapid pipeline done"
