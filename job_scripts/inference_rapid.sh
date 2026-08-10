#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           inference_rapid.sh
###############################################################################

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --time=1-00:00:00
#SBATCH --partition=qcb
#SBATCH --account=jazlynmo_738
#SBATCH --nodes=1
#SBATCH --output=/home1/karatas/logs/rapid/%A_%a.%x.out
#SBATCH --error=/home1/karatas/logs/rapid/%A_%a.%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=karatas@usc.edu

set -euo pipefail
ml purge
ml gcc/13.3.0 htslib/1.19.1
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh

: "${SLURM_ARRAY_TASK_ID:?ERROR: run as a Slurm array}"

# input
min_len="$1"
window_size="$2"
num_runs="$3"
num_successes="$4"
vcf_dir="$5"
rapid_map_dir="$6"
outdir="$7"
num_reps="$8"
chr="$9"
shift 9

# read TMAPs until "--", then read imaps until "--", thenread pops
tmaps=()
while [[ $# -gt 0 && "$1" != "--" ]]; do
        tmaps+=( "$1" )
        shift
done
[[ $# -gt 0 && "$1" == "--" ]] || { echo "ERROR: expected '--' after tmaps" >&2; exit 1; }
shift 1

imaps=()
while [[ $# -gt 0 && "$1" != "--" ]]; do
        imaps+=( "$1" )
        shift
done
[[ $# -gt 0 && "$1" == "--" ]] || { echo "ERROR: expected '--' after imaps" >&2; exit 1; }
shift 1

pops=( "$@" )

[[ ${#tmaps[@]} -gt 0 ]] || { echo "ERROR: no tmaps provided" >&2; exit 1; }
[[ ${#imaps[@]} -gt 0 ]] || { echo "ERROR: no imaps provided" >&2; exit 1; }
[[ ${#pops[@]} -gt 0 ]] || { echo "ERROR: no pops provided" >&2; exit 1; }

# decode task id to which tmap + imap + replicate
tid="${SLURM_ARRAY_TASK_ID}"
tid0=$(( tid - 1 ))

rep=$(( tid0 % num_reps + 1 ))
tmp_id=$(( tid0 / num_reps ))

imap_idx=$(( tmp_id % ${#imaps[@]} ))
tmap_idx=$(( tmp_id / ${#imaps[@]} ))

if (( tmap_idx < 0 || tmap_idx >= ${#tmaps[@]} )); then
        echo "ERROR: task ${tid} produced tmap_idx=${tmap_idx}, but n_maps=${#tmaps[@]}" >&2
        exit 1
fi
if (( imap_idx < 0 || imap_idx >= ${#imaps[@]} )); then
        echo "ERROR: task ${tid} produced imap_idx=${imap_idx}, but n_maps=${#imaps[@]}" >&2
        exit 1
fi

tmap="${tmaps[$tmap_idx]}"
imap="${imaps[$imap_idx]}"

# for each pop, infer ibd
ibd_dir="${outdir}/rapid_sim_tmap${tmap}_imap${imap}_d${min_len}_w${window_size}_r${num_runs}_s${num_successes}"
mkdir -p "$ibd_dir"
log_msg "ibd_dir defined as ${ibd_dir}"
for pop in "${pops[@]}"; do
        vcf="${vcf_dir}/${tmap}_${rep}_chr${chr}_${pop}.vcf.gz"
        [[ -s "$vcf" ]] || { log_msg "vcf not found: $vcf"; exit 1; }

        genetic_map="${rapid_map_dir}/pop${pop}_chr${chr}_rep${rep}_tmap${tmap}_imap${imap}.map"
        [[ -s "${genetic_map}" ]] || {
                log_msg "cached rapid map missing: ${genetic_map}"
                exit 1
        }

        ibd_out="${ibd_dir}/${pop}_chr${chr}_rep${rep}_out"
        ibd_bed="${ibd_dir}/${pop}_chr${chr}_rep${rep}_ibd.bed.gz"
        if [[ -s "$ibd_bed" ]]; then
                log_msg "${ibd_bed} already found for ${pop}; continue"
                continue
        fi

        # call ibd with rapid
        log_msg "Calling IBD of $pop for rep${rep} with rapid using...
        min_length_in_cM = ${min_len}
        window_size = ${window_size}
        num_runs = ${num_runs}
        num_successes = ${num_successes}
        input_file_vcf_compressed = ${vcf}
        genetic_mapping_file = ${genetic_map}
        output_folder = $ibd_out"

        rm -rf "$ibd_out"
        "${HOME}/software/RaPID/RaPID_v.1.7" \
                -i "$vcf" \
                -g "$genetic_map" \
                -d "$min_len" \
                -w "$window_size" \
                -r "$num_runs" \
                -s "$num_successes" \
                -o "$ibd_out"

        # exit nonzero if rapid did not produce an output
        if [[ ! -s "${ibd_out}/results.max.gz" ]]; then
                log_msg "no output for $pop rep$rep"
                exit 1
        fi

        # early exit if no ibd segments were called
        if [[ $(zcat "${ibd_out}/results.max.gz" | wc -l) -eq 0 ]]; then
                log_msg "no ibd segments detected for $pop on rep$rep; continue"
                continue
        fi

        # convert ibd calls to bed
        tmpdir="${SCRATCH:-/tmp}/bedtmp_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID:-$$}"
        mkdir -p "$tmpdir"
        tmp_bed="${ibd_bed}.tmp"
        zcat "${ibd_out}/results.max.gz" | \
        awk -F'\t' 'BEGIN{OFS="\t"}   
                $6 ~ /^-?[0-9]+$/ && $7 ~ /^-?[0-9]+$/ && $6 >= 0 && $7 >= $6 {
                print $1,$6,$7,$2,$4,$3,$5,$8,$9,$10
                }' | \
        sort -k1,1 -k2,2n -k3,3n -T "$tmpdir" -S 2G | \
        gzip > "$tmp_bed"

        mv -f "$tmp_bed" "$ibd_bed"
        rm -rf "$tmpdir"
done

log_msg "called ibd for tmap=${tmap}, imap=${imap}, rep=${rep}"
