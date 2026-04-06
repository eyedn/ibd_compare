#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           filter_format_inference.sh
###############################################################################

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=1-00:00:00
#SBATCH --partition=qcb
#SBATCH --account=jazlynmo_738
#SBATCH --nodes=1
#SBATCH --output=/home1/karatas/logs/infFF.%A_%a.out
#SBATCH --error=/home1/karatas/logs/infFF.%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=karatas@usc.edu

set -euo pipefail
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh

: "${SLURM_ARRAY_TASK_ID:?ERROR: run as a Slurm array}"

# hapibd params (to find correct ibd file)
MIN_SEED=2.0
MAX_GAP=1000
MIN_EXTEND=1.0
MIN_OUTPUT=2.0
MIN_MARKERS=100
MIN_MAC=2

# rapid params (to find correct ibd file)
MIN_LEN=2
WINDOW_SIZE=300
NUM_RUNS=10
NUM_SUCCESSES=2

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
shift 10

# read TMAPs until "--", then read IMAPS until "--", thenread POPS
TMAPS=()
while [[ $# -gt 0 && "$1" != "--" ]]; do
        TMAPS+=( "$1" )
        shift
done
[[ $# -gt 0 && "$1" == "--" ]] || { echo "ERROR: expected '--' after TMAPS" >&2; exit 1; }
shift 1

IMAPS=()
while [[ $# -gt 0 && "$1" != "--" ]]; do
        IMAPS+=( "$1" )
        shift
done
[[ $# -gt 0 && "$1" == "--" ]] || { echo "ERROR: expected '--' after IMAPS" >&2; exit 1; }
shift 1

POPS=( "$@" )

[[ ${#TMAPS[@]} -gt 0 ]] || { echo "ERROR: no TMAPS provided" >&2; exit 1; }
[[ ${#IMAPS[@]} -gt 0 ]] || { echo "ERROR: no IMAPS provided" >&2; exit 1; }
[[ ${#POPS[@]} -gt 0 ]] || { echo "ERROR: no POPS provided" >&2; exit 1; }

# decode task id to which tmap + imap + replicate/seed
tid="${SLURM_ARRAY_TASK_ID}"
tid0=$(( tid - 1 ))

seed_id=$(( tid0 % NUM_REPS + 1 ))
tmp_id=$(( tid0 / NUM_REPS ))

imap_idx=$(( tmp_id % ${#IMAPS[@]} ))
tmap_idx=$(( tmp_id / ${#IMAPS[@]} ))

if (( tmap_idx < 0 || tmap_idx >= ${#TMAPS[@]} )); then
        echo "ERROR: task ${tid} produced tmap_idx=${tmap_idx}, but n_maps=${#TMAPS[@]}" >&2
        exit 1
fi
if (( imap_idx < 0 || imap_idx >= ${#IMAPS[@]} )); then
        echo "ERROR: task ${tid} produced imap_idx=${imap_idx}, but n_maps=${#IMAPS[@]}" >&2
        exit 1
fi

tmap="${TMAPS[$tmap_idx]}"
imap="${IMAPS[$imap_idx]}"
seed="${seed_id}"

mkdir -p "$COMMON_FF_DIR"

# loop through each pop
for pop in "${POPS[@]}"; do
    vcf="${VCF_DIR}/${tmap}_${seed}_${pop}.vcf.gz"

    # create standardized tmap label for final ibd file
    if [[ "$tmap" == *HapMapII_GRCh38* ]]; then 
        tmap_label="hapmap"
    elif [[ "$tmap" == *DeCodeSexAveraged_GRCh38* ]]; then
        tmap_label="decode"
    elif [[ "$tmap" == *DefaultUniform* ]]; then
        tmap_label="flat"
    fi

    # define paths for the true gmaps for cM interpolation of ibd str/end;
    # create standardized imap label for final ibd file
    if [[ "$imap" == *HapMapII_GRCh38* ]]; then 
        imap_label="hapmap"
        hapibd_gmap="${REF_DIR}/HapMapII_GRCh38/plink.chr${CHR}.GRCh38.UCSC.format.map"
    elif [[ "$imap" == *DeCodeSexAveraged_GRCh38* ]]; then
        imap_label="decode"
        hapibd_gmap="${REF_DIR}/deCODE/addendum/sex-averaged_gmaps/sex_av_chr${CHR}_hg38.gmap"
    elif [[ "$imap" == *Flat1Mb1cM* ]]; then
        imap_label="1Mb1cM"
        hapibd_gmap="${REF_DIR}/GRCh38.p14/recomb_maps_1Mb1cM/chr${CHR}.map"
    fi
    rapid_gmap_dir="${OUTDIR}/rapid_sim_tmap${tmap}_imap${imap}_d${MIN_LEN}_w${WINDOW_SIZE}_r${NUM_RUNS}_s${NUM_SUCCESSES}/maps"
    rapid_gmap="${rapid_gmap_dir}/chr${CHR}.${pop}.${seed}.RaPIDfinal.map"

    # general paths
    hapibd_dir_base="hapibd_sim_tmap${tmap}_imap${imap}_seed${MIN_SEED}_gap${MAX_GAP}_extend${MIN_EXTEND}_out${MIN_OUTPUT}_mark${MIN_MARKERS}_mac${MIN_MAC}"
    hapibd_dir="${OUTDIR}/${hapibd_dir_base}"

    rapid_dir_base="rapid_sim_tmap${tmap}_imap${imap}_d${MIN_LEN}_w${WINDOW_SIZE}_r${NUM_RUNS}_s${NUM_SUCCESSES}"
    rapid_dir="${OUTDIR}/${rapid_dir_base}"

    # paths related to ibd-ends
    hapibd_nonbed="${hapibd_dir}/${pop}_chr${CHR}_seed${seed}.ibd.gz"
    hapibd_ibdends_prefix="${hapibd_dir}/${pop}_chr${CHR}_seed${seed}.ibdends"

    rapid_nonbed_dir="${rapid_dir}/${pop}_chr${CHR}_seed${seed}_out"
    rapid_nonbed="${rapid_nonbed_dir}/results.max.gz"
    rapid_nonbed_hapibdformat="${rapid_nonbed_dir}/results.max.hapibdformat.gz"
    rapid_ibdends_prefix="${rapid_dir}/${pop}_chr${CHR}_seed${seed}.ibdends"

    # inputs to filter/formatting by R script
    hapibd_ibd="${hapibd_dir}/${pop}_chr${CHR}_seed${seed}_ibd.bed.gz"
    hapibd_ibdends_ibd="${hapibd_ibdends_prefix}_ibd.ibdends.bed.gz"

    rapid_ibd="${rapid_dir}/${pop}_chr${CHR}_seed${seed}_ibd.bed.gz"
    rapid_ibdends_ibd="${rapid_ibdends_prefix}_ibd.ibdends.bed.gz"

    # outputs from filter/formatting R script
    hapibd_ff_dir="${COMMON_FF_DIR}/${hapibd_dir_base}"
    hapibd_out="${hapibd_ff_dir}/hapibd_${tmap_label}_${imap_label}_${pop}_${seed}.bed"
    hapibd_hot_out="${hapibd_ff_dir}/hapibd_${tmap_label}_${imap_label}_${pop}_${seed}.hot.bed"
    hapibd_ibdends_out="${hapibd_ff_dir}/hapibd_${tmap_label}_${imap_label}_${pop}_${seed}.ibdends.bed"
    hapibd_ibdendshot_out="${hapibd_ff_dir}/hapibd_${tmap_label}_${imap_label}_${pop}_${seed}.ibdendshot.bed"
    
    rapid_ff_dir="${COMMON_FF_DIR}/${rapid_dir_base}"
    rapid_out="${rapid_ff_dir}/rapid_${tmap_label}_${imap_label}_${pop}_${seed}.bed"
    rapid_hot_out="${rapid_ff_dir}/rapid_${tmap_label}_${imap_label}_${pop}_${seed}.hot.bed"
    rapid_ibdends_out="${rapid_ff_dir}/rapid_${tmap_label}_${imap_label}_${pop}_${seed}.ibdends.bed"
    rapid_ibdendshot_out="${rapid_ff_dir}/rapid_${tmap_label}_${imap_label}_${pop}_${seed}.ibdendshot.bed"
    
    [[ -s "$vcf" ]] || { echo "ERROR: missing vcf: $vcf" >&2; exit 1; }
    [[ -s "$hapibd_nonbed" ]] || { echo "ERROR: missing hapibd ibd.gz: $hapibd_nonbed" >&2; exit 1; }
    [[ -s "$rapid_nonbed" ]] || { echo "ERROR: missing rapid results: $rapid_nonbed" >&2; exit 1; }

    mkdir -p "$hapibd_ff_dir"
    mkdir -p "$rapid_ff_dir"
    
    # filter hapibd calls by ibd-ends
    module purge
    ml legacy/CentOS7 jdk/17.0.5
    log_msg "filtering hapibd segments using ibd-ends"
    java -Xmx46G -jar ~/software/ibd-ends/ibd-ends.jar \
        gt="$vcf" \
        ibd="$hapibd_nonbed" \
        map="$hapibd_gmap" \
        nthreads="${SLURM_CPUS_PER_TASK:-1}" \
        out="$hapibd_ibdends_prefix"
    log_msg "generating ${hapibd_ibdends_ibd} from ${hapibd_ibdends_prefix}.ibd.gz"
    zcat "${hapibd_ibdends_prefix}.ibd.gz" | \
        awk 'BEGIN{OFS="\t"} NR>1 {print $5,$10,$11,$1,$2,$3,$4,$12}' | \
        gzip -c > "$hapibd_ibdends_ibd"
    
    # filter rapid calls by ibd-ends (note, non-rapid format map needed here)
    log_msg "filtering rapid segmends using ibd-ends"
    log_msg "generating ${rapid_nonbed_hapibdformat} (hapibd format) from ${rapid_nonbed} (rapid format)"
    zcat "$rapid_nonbed" | \
        awk 'BEGIN{OFS="\t"}
            ($6 >= 0 && $7 >= 0 && $6 <= $7) {
                print $2, $4+1, $3, $5+1, $1, $6, $7, $8
            }' | \
        gzip -c > "$rapid_nonbed_hapibdformat"
    java -Xmx46G -jar ~/software/ibd-ends/ibd-ends.jar \
        gt="$vcf" \
        ibd="${rapid_nonbed_hapibdformat}" \
        map="$hapibd_gmap" \
        nthreads="${SLURM_CPUS_PER_TASK:-1}" \
        out="$rapid_ibdends_prefix"
    log_msg "generating ${rapid_ibdends_ibd} (hapibd format) from ${rapid_ibdends_prefix}.ibd.gz (ibdends format)"
    zcat "${rapid_ibdends_prefix}.ibd.gz" | \
        awk 'BEGIN{OFS="\t"} NR>1 {print $5,$10,$11,$1,$2,$3,$4,$12}' | \
        gzip -c > "$rapid_ibdends_ibd"

    # format hapibd calls and filter out hot bins (bins with large no. ibd)
    module purge
    ml rstats/4.5.1
    log_msg "format and filtering hapibd calls"
    Rscript "/home1/karatas/ibd_compare/job_scripts/filter_and_format.R" \
        "$MIN_BP_LEN" \
        "$MIN_CM_LEN" \
        "$HOT_BIN_SIZE" \
        "$HOT_BIN_SD" \
        "$hapibd_ibd" \
        "$hapibd_gmap" \
        "$hapibd_out" \
        "$hapibd_hot_out" \
        "TRUE"
    gzip -f "$hapibd_out"
    gzip -f "$hapibd_hot_out"

    # format ibdends-hapibd calls and filter out hot bins (bins with large no. ibd)
    log_msg "format and filtering hapibd calls adjusted by ibdends"
    Rscript "/home1/karatas/ibd_compare/job_scripts/filter_and_format.R" \
        "$MIN_BP_LEN" \
        "$MIN_CM_LEN" \
        "$HOT_BIN_SIZE" \
        "$HOT_BIN_SD" \
        "$hapibd_ibdends_ibd" \
        "$hapibd_gmap" \
        "$hapibd_ibdends_out" \
        "$hapibd_ibdendshot_out" \
        "TRUE"
    gzip -f "$hapibd_ibdends_out"
    gzip -f "$hapibd_ibdendshot_out"

    # format rapid calls and filter out hot bins (bins with large no. ibd)
    log_msg "format and filtering rapid calls"
    Rscript "/home1/karatas/ibd_compare/job_scripts/filter_and_format.R" \
        "$MIN_BP_LEN" \
        "$MIN_CM_LEN" \
        "$HOT_BIN_SIZE" \
        "$HOT_BIN_SD" \
        "$rapid_ibd" \
        "$rapid_gmap" \
        "$rapid_out" \
        "$rapid_hot_out" \
        "FALSE"
    gzip -f "$rapid_out" 
    gzip -f "$rapid_hot_out"

    # format ibdends-rapid calls and filter out hot bins (bins with large no. ibd)
    log_msg "format and filtering rapid calls adjusted by ibdends (note, hapibd bed format here with non-rapid gmap)"
    Rscript "/home1/karatas/ibd_compare/job_scripts/filter_and_format.R" \
        "$MIN_BP_LEN" \
        "$MIN_CM_LEN" \
        "$HOT_BIN_SIZE" \
        "$HOT_BIN_SD" \
        "$rapid_ibdends_ibd" \
        "$hapibd_gmap" \
        "$rapid_ibdends_out" \
        "$rapid_ibdendshot_out" \
        "TRUE"
    gzip -f "$rapid_ibdends_out"
    gzip -f "$rapid_ibdendshot_out"
done

log_msg "done"
