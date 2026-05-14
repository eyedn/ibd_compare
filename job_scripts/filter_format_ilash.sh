#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           filter_format_ilash.sh
###############################################################################

#SBATCH --time=1-00:00:00
#SBATCH --partition=qcb
#SBATCH --account=jazlynmo_738
#SBATCH --nodes=1
#SBATCH --output=/home1/karatas/logs/ilaFF.%A_%a.out
#SBATCH --error=/home1/karatas/logs/ilaFF.%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=karatas@usc.edu

set -euo pipefail
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh

: "${IBDENDS_XMX:?ERROR: IBDENDS_XMX must be set by submit wrapper}"
: "${IBDENDS_THREADS:?ERROR: IBDENDS_THREADS must be set by submit wrapper}"
: "${IBDENDS_MIN_MAF:?ERROR: IBDENDS_MIN_MAF must be set by submit wrapper}"

ibdends_xmx="${IBDENDS_XMX}"
ibdends_threads="${IBDENDS_THREADS}"
ibdends_min_maf="${IBDENDS_MIN_MAF}"

: "${SLURM_ARRAY_TASK_ID:?ERROR: run as a Slurm array}"

# inputs
outdir="$1"
vcf_dir="$2"
ilash_map_dir="$3"
common_ff_dir="$4"
chr="$5"
min_bp_len="$6"
min_cm_len="$7"
hot_bin_size="$8"
hot_bin_sd="$9"
num_reps="${10}"
ilash_auto_slice="${11}"
ilash_min_length="${12}"
ilash_slice_length="${13}"
ilash_cm_overlap="${14}"
shift 14

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
if (( tmap_idx < 0 || tmap_idx >= ${#tmaps[@]} )); then
    echo "ERROR: task ${tid} produced tmap_idx=${tmap_idx}" >&2
    exit 1
fi
if (( imap_idx < 0 || imap_idx >= ${#imaps[@]} )); then
    echo "ERROR: task ${tid} produced imap_idx=${imap_idx}" >&2
    exit 1
fi

tmap="${tmaps[$tmap_idx]}"
imap="${imaps[$imap_idx]}"
seed="${seed_id}"
pop="${pops[$pop_idx]}"

mkdir -p "${common_ff_dir}"

if [[ "$tmap" == *HapMapII_GRCh38* ]]; then tmap_label="hapmap"
elif [[ "$tmap" == *DeCodeSexAveraged_GRCh38* ]]; then tmap_label="decode"
elif [[ "$tmap" == *DefaultUniform* ]]; then tmap_label="flat"
else tmap_label="$tmap"
fi

if [[ "$imap" == *HapMapII_GRCh38* ]]; then
    imap_label="hapmap"
elif [[ "$imap" == *DeCodeSexAveraged_GRCh38* ]]; then
    imap_label="decode"
elif [[ "$imap" == *Flat1Mb1cM* ]]; then
    imap_label="1Mb1cM"
else
    imap_label="$imap"
fi

ilash_dir_base="ilash_sim_tmap${tmap}_imap${imap}_auto${ilash_auto_slice}_min${ilash_min_length}_slice${ilash_slice_length}_overlap${ilash_cm_overlap}"
ilash_dir="${outdir}/${ilash_dir_base}"

vcf="${vcf_dir}/${tmap}_${seed}_${pop}.vcf.gz"
ilash_gmap="${ilash_map_dir}/pop${pop}_chr${chr}_seed${seed}_tmap${tmap}_imap${imap}.map"
ilash_nonbed="${ilash_dir}/${pop}_chr${chr}_seed${seed}.ibd.gz"
ilash_ibd="${ilash_dir}/${pop}_chr${chr}_seed${seed}_ibd.bed.gz"
ilash_ibdends_prefix="${ilash_dir}/${pop}_chr${chr}_seed${seed}.ibdends"
ilash_ibdends_ibd="${ilash_ibdends_prefix}_ibd.ibdends.bed.gz"
ilash_hot_hapibdformat="${ilash_dir}/${pop}_chr${chr}_seed${seed}.hot.hapibdformat.gz"
ilash_hot_ibdends_prefix="${ilash_dir}/${pop}_chr${chr}_seed${seed}.hot.ibdends"
ilash_hot_ibdends_ibd="${ilash_hot_ibdends_prefix}_ibd.ibdends.bed.gz"

ilash_ff_dir="${common_ff_dir}/${ilash_dir_base}"
ilash_out="${ilash_ff_dir}/ilash_${tmap_label}_${imap_label}_${pop}_${seed}.bed"
ilash_hot_out="${ilash_ff_dir}/ilash_${tmap_label}_${imap_label}_${pop}_${seed}.hot.bed"
ilash_ibdends_out="${ilash_ff_dir}/ilash_${tmap_label}_${imap_label}_${pop}_${seed}.ibdends.bed"
ilash_hotibdends_out="${ilash_ff_dir}/ilash_${tmap_label}_${imap_label}_${pop}_${seed}.hotibdends.bed"
ilash_ibdendshot_out="${ilash_ff_dir}/ilash_${tmap_label}_${imap_label}_${pop}_${seed}.ibdendshot.bed"

[[ -s "${vcf}" ]] || { echo "ERROR: missing vcf: ${vcf}" >&2; exit 1; }
[[ -s "${ilash_gmap}" ]] || { echo "ERROR: missing ilash map: ${ilash_gmap}" >&2; exit 1; }
[[ -s "${ilash_nonbed}" ]] || { echo "ERROR: missing ilash ibd.gz: ${ilash_nonbed}" >&2; exit 1; }
[[ -s "${ilash_ibd}" ]] || { echo "ERROR: missing ilash bed.gz: ${ilash_ibd}" >&2; exit 1; }

mkdir -p "${ilash_ff_dir}"

module purge
ml openjdk/25.0.2
log_msg "ilash: ibd-ends runtime
ibdends_xmx=${ibdends_xmx}
ibdends_threads=${ibdends_threads}"
log_msg "ilash: filtering segments using ibd-ends (pop=${pop} seed=${seed} tmap=${tmap_label} imap=${imap_label})"
java -Xmx"${ibdends_xmx}" -jar ~/software/ibd-ends/ibd-ends.jar \
    gt="${vcf}" \
    ibd="${ilash_nonbed}" \
    map="${ilash_gmap}" \
    min-maf="${ibdends_min_maf}" \
    nthreads="${ibdends_threads}" \
    out="${ilash_ibdends_prefix}"

log_msg "ilash: generating ${ilash_ibdends_ibd} from ${ilash_ibdends_prefix}.ibd.gz"
zcat "${ilash_ibdends_prefix}.ibd.gz" | \
    awk 'BEGIN{OFS="\t"} NR>1 {print $5,$10,$11,$1,$2,$3,$4,$12}' | \
    gzip -c > "${ilash_ibdends_ibd}"

module purge
ml rstats/4.5.1

log_msg "ilash: format+filter raw calls"
Rscript "/home1/karatas/ibd_compare/job_scripts/filter_and_format.R" \
    "${min_bp_len}" "${min_cm_len}" "${hot_bin_size}" "${hot_bin_sd}" \
    "${ilash_ibd}" "${ilash_gmap}" \
    "${ilash_out}" "${ilash_hot_out}" \
    "TRUE" "TRUE"

log_msg "ilash: converting hot-filtered calls to hapibd format"
Rscript "/home1/karatas/ibd_compare/job_scripts/convert_std_bed_to_hapibd.R" \
    "${ilash_hot_out}" "${ilash_hot_hapibdformat}"

module purge
ml openjdk/25.0.2

log_msg "ilash: filtering hot-filtered segments using ibd-ends"
java -Xmx"${ibdends_xmx}" -jar ~/software/ibd-ends/ibd-ends.jar \
    gt="${vcf}" \
    ibd="${ilash_hot_hapibdformat}" \
    map="${ilash_gmap}" \
    min-maf="${ibdends_min_maf}" \
    nthreads="${ibdends_threads}" \
    out="${ilash_hot_ibdends_prefix}"

log_msg "ilash: generating ${ilash_hot_ibdends_ibd} from ${ilash_hot_ibdends_prefix}.ibd.gz"
zcat "${ilash_hot_ibdends_prefix}.ibd.gz" | \
    awk 'BEGIN{OFS="\t"} NR>1 {print $5,$10,$11,$1,$2,$3,$4,$12}' | \
    gzip -c > "${ilash_hot_ibdends_ibd}"

module purge
ml rstats/4.5.1

log_msg "ilash: format+filter ibd-ends adjusted calls"
Rscript "/home1/karatas/ibd_compare/job_scripts/filter_and_format.R" \
    "${min_bp_len}" "${min_cm_len}" "${hot_bin_size}" "${hot_bin_sd}" \
    "${ilash_ibdends_ibd}" "${ilash_gmap}" \
    "${ilash_ibdends_out}" "${ilash_ibdendshot_out}" \
    "TRUE" "TRUE"

log_msg "ilash: formatting hot-filtered ibd-ends adjusted calls"
Rscript "/home1/karatas/ibd_compare/job_scripts/filter_and_format.R" \
    "${min_bp_len}" "${min_cm_len}" "${hot_bin_size}" "${hot_bin_sd}" \
    "${ilash_hot_ibdends_ibd}" "${ilash_gmap}" \
    "${ilash_hotibdends_out}" "/tmp/ilash_hotibdends_unused.bed" \
    "TRUE" "FALSE"

gzip -f \
    "${ilash_out}" \
    "${ilash_hot_out}" \
    "${ilash_ibdends_out}" \
    "${ilash_hotibdends_out}" \
    "${ilash_ibdendshot_out}"

log_msg "ilash pipeline done"
