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

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16G
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
outdir="$1"
vcf_dir="$2"
common_ff_dir="$3"
ref_dir="$4"
chr="$5"
min_bp_len="$6"
min_cm_len="$7"
hot_bin_size="$8"
hot_bin_sd="$9"
num_reps="${10}"

# rapid params (from submit script)
rap_min_len="${11}"
rap_window_size="${12}"
rap_num_runs="${13}"
rap_num_successes="${14}"
shift 14

# read TMAPs until "--", then imaps until "--", then pops
tmaps=()
while [[ $# -gt 0 && "$1" != "--" ]]; do tmaps+=( "$1" ); shift; done
[[ $# -gt 0 && "$1" == "--" ]] || { echo "ERROR: expected '--' after tmaps" >&2; exit 1; }
shift

imaps=()
while [[ $# -gt 0 && "$1" != "--" ]]; do imaps+=( "$1" ); shift; done
[[ $# -gt 0 && "$1" == "--" ]] || { echo "ERROR: expected '--' after imaps" >&2; exit 1; }
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
    echo "ERROR: task ${tid} produced pop_idx=${pop_idx}, n_pops=${#pops[@]}" >&2
    exit 1
fi
if (( tmap_idx < 0 || tmap_idx >= ${#tmaps[@]} )); then
    echo "ERROR: task ${tid} produced tmap_idx=${tmap_idx}, n_tmaps=${#tmaps[@]}" >&2
    exit 1
fi
if (( imap_idx < 0 || imap_idx >= ${#imaps[@]} )); then
    echo "ERROR: task ${tid} produced imap_idx=${imap_idx}, n_imaps=${#imaps[@]}" >&2
    exit 1
fi

tmap="${tmaps[$tmap_idx]}"
imap="${imaps[$imap_idx]}"
seed="${seed_id}"
pop="${pops[$pop_idx]}"

mkdir -p "$common_ff_dir"

# standardized labels
if [[ "$tmap" == *HapMapII_GRCh38* ]]; then tmap_label="hapmap"
elif [[ "$tmap" == *DeCodeSexAveraged_GRCh38* ]]; then tmap_label="decode"
elif [[ "$tmap" == *DefaultUniform* ]]; then tmap_label="flat"
else tmap_label="$tmap"
fi

# NOTE: ibd-ends needs a NON-RaPID-style map; you were using hapibd_gmap for this
if [[ "$imap" == *HapMapII_GRCh38* ]]; then
    imap_label="hapmap"
    hapibd_gmap="${ref_dir}/HapMapII_GRCh38/plink.chr${chr}.GRCh38.UCSC.format.map"
elif [[ "$imap" == *DeCodeSexAveraged_GRCh38* ]]; then
    imap_label="decode"
    hapibd_gmap="${ref_dir}/deCODE/addendum/sex-averaged_gmaps/sex_av_chr${chr}_hg38.gmap"
elif [[ "$imap" == *Flat1Mb1cM* ]]; then
    imap_label="1Mb1cM"
    hapibd_gmap="${ref_dir}/GRCh38.p14/recomb_maps_1Mb1cM/chr${chr}.map"
else
    imap_label="$imap"
    hapibd_gmap=""
fi
[[ -n "$hapibd_gmap" && -s "$hapibd_gmap" ]] || { echo "ERROR: missing hapibd_gmap: $hapibd_gmap" >&2; exit 1; }

# rapid dir base (use const params)
rapid_dir_base="rapid_sim_tmap${tmap}_imap${imap}_d${rap_min_len}_w${rap_window_size}_r${rap_num_runs}_s${rap_num_successes}"
rapid_dir="${outdir}/${rapid_dir_base}"

vcf="${vcf_dir}/${tmap}_${seed}_${pop}.vcf.gz"

rapid_nonbed_dir="${rapid_dir}/${pop}_chr${chr}_seed${seed}_out"
rapid_nonbed="${rapid_nonbed_dir}/results.max.gz"
rapid_nonbed_hapibdformat="${rapid_nonbed_dir}/results.max.hapibdformat.gz"
rapid_ibdends_prefix="${rapid_dir}/${pop}_chr${chr}_seed${seed}.ibdends"
rapid_ibdends_ibd="${rapid_ibdends_prefix}_ibd.ibdends.bed.gz"
rapid_hot_hapibdformat="${rapid_nonbed_dir}/results.max.hot.hapibdformat.gz"
rapid_hot_ibdends_prefix="${rapid_dir}/${pop}_chr${chr}_seed${seed}.hot.ibdends"
rapid_hot_ibdends_ibd="${rapid_hot_ibdends_prefix}_ibd.ibdends.bed.gz"

# bed made from previous step
rapid_ibd="${rapid_dir}/${pop}_chr${chr}_seed${seed}_ibd.bed.gz"

# rapid map for cM interpolation during formatting
rapid_gmap_dir="${rapid_dir}/maps"
rapid_gmap="${rapid_gmap_dir}/chr${chr}.${pop}.${seed}.RaPIDfinal.map"

rapid_ff_dir="${common_ff_dir}/${rapid_dir_base}"
rapid_out="${rapid_ff_dir}/rapid_${tmap_label}_${imap_label}_${pop}_${seed}.bed"
rapid_hot_out="${rapid_ff_dir}/rapid_${tmap_label}_${imap_label}_${pop}_${seed}.hot.bed"
rapid_ibdends_out="${rapid_ff_dir}/rapid_${tmap_label}_${imap_label}_${pop}_${seed}.ibdends.bed"
rapid_hotibdends_out="${rapid_ff_dir}/rapid_${tmap_label}_${imap_label}_${pop}_${seed}.hotibdends.bed"
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
java -Xmx48G -jar ~/software/ibd-ends/ibd-ends.jar \
    gt="$vcf" \
    ibd="$rapid_nonbed_hapibdformat" \
    map="$hapibd_gmap" \
    nthreads="4" \
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
    "$min_bp_len" "$min_cm_len" "$hot_bin_size" "$hot_bin_sd" \
    "$rapid_ibd" "$rapid_gmap" \
    "$rapid_out" "$rapid_hot_out" \
    "FALSE" "TRUE"

log_msg "rapid: converting hot-filtered calls to hapibd format"
Rscript "/home1/karatas/ibd_compare/job_scripts/convert_std_bed_to_hapibd.R" \
    "$rapid_hot_out" "$rapid_hot_hapibdformat"

log_msg "rapid: format+filter ibd-ends adjusted calls (hapibd-style map here)"
Rscript "/home1/karatas/ibd_compare/job_scripts/filter_and_format.R" \
    "$min_bp_len" "$min_cm_len" "$hot_bin_size" "$hot_bin_sd" \
    "$rapid_ibdends_ibd" "$hapibd_gmap" \
    "$rapid_ibdends_out" "$rapid_ibdendshot_out" \
    "TRUE" "TRUE"

module purge
ml legacy/CentOS7 jdk/17.0.5

log_msg "rapid: filtering hot-filtered segments using ibd-ends"
java -Xmx48G -jar ~/software/ibd-ends/ibd-ends.jar \
    gt="$vcf" \
    ibd="$rapid_hot_hapibdformat" \
    map="$hapibd_gmap" \
    nthreads="4" \
    out="$rapid_hot_ibdends_prefix"

log_msg "rapid: generating ${rapid_hot_ibdends_ibd} from ${rapid_hot_ibdends_prefix}.ibd.gz"
zcat "${rapid_hot_ibdends_prefix}.ibd.gz" | \
    awk 'BEGIN{OFS="\t"} NR>1 {print $5,$10,$11,$1,$2,$3,$4,$12}' | \
    gzip -c > "$rapid_hot_ibdends_ibd"

module purge
ml rstats/4.5.1

log_msg "rapid: formatting hot-filtered ibd-ends adjusted calls"
Rscript "/home1/karatas/ibd_compare/job_scripts/filter_and_format.R" \
    "$min_bp_len" "$min_cm_len" "$hot_bin_size" "$hot_bin_sd" \
    "$rapid_hot_ibdends_ibd" "$hapibd_gmap" \
    "$rapid_hotibdends_out" "/tmp/rapid_hotibdends_unused.bed" \
    "TRUE" "FALSE"

gzip -f \
    "$rapid_out" \
    "$rapid_hot_out" \
    "$rapid_ibdends_out" \
    "$rapid_hotibdends_out" \
    "$rapid_ibdendshot_out"

log_msg "rapid pipeline done"
