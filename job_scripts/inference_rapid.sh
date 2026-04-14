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
#SBATCH --output=/home1/karatas/logs/ibd_array.%A_%a.out
#SBATCH --error=/home1/karatas/logs/ibd_array.%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=karatas@usc.edu

set -euo pipefail
ml purge
ml gcc/13.3.0 htslib/1.19.1 bcftools/1.19 conda
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh
source /apps/conda/miniforge3/25.3.0/etc/profile.d/conda.sh
conda activate rpgg_proj
export PATH="/home1/karatas/.conda/envs/rpgg_proj/bin:$PATH"
log_msg "using python: $(which python)"
python -c "import sys; print(sys.executable)"
python -c "import numpy; print('NumPy OK:', numpy.__version__)"

: "${SLURM_ARRAY_TASK_ID:?ERROR: run as a Slurm array}"

# input
min_len="$1"
window_size="$2"
num_runs="$3"
num_successes="$4"
outdir="$5"
vcf_dir="$6"
num_reps="$7"
chr="$8"
shift 8

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

# decode task id to which tmap + imap + replicate/seed
tid="${SLURM_ARRAY_TASK_ID}"
tid0=$(( tid - 1 ))

seed_id=$(( tid0 % num_reps + 1 ))
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
seed="${seed_id}"

# define intermediate (chr, bp, cM) that will be filtered and interpolated
if [[ "$imap" == "HapMapII_GRCh38" ]]; then
        int_genetic_map="/home1/karatas/references/RaPID_format_maps/HapMapII_GRCh38/plink.chr${chr}.GRCh38.UCSC.reformated.map"
elif [[ "$imap" == "DeCodeSexAveraged_GRCh38" ]]; then
        int_genetic_map="/home1/karatas/references/RaPID_format_maps/deCODE/chr${chr}.reformated.map"
elif [[ "$imap" == "Flat1Mb1cM" ]]; then
        int_genetic_map="/home1/karatas/references/RaPID_format_maps/1Mb1cM/chr${chr}.reformated.map"
else
        echo "ERROR: unknown imap=${imap}" >&2
        exit 1
fi

# for each pop, infer ibd
ibd_dir="${outdir}/rapid_sim_tmap${tmap}_imap${imap}_d${min_len}_w${window_size}_r${num_runs}_s${num_successes}"
mkdir -p "$ibd_dir"
mkdir -p "$ibd_dir/maps"
log_msg "ibd_dir defined as ${ibd_dir}"
for pop in "${pops[@]}"; do
        vcf="${vcf_dir}/${tmap}_${seed}_${pop}.vcf.gz"
        [[ -s "$vcf" ]] || { log_msg "vcf not found: $vcf"; exit 1; }

        # define filtered map and final genetic map paths after rapid filtering
        # and interpolating for this run
        filt_genetic_map="${ibd_dir}/maps/chr${chr}.${pop}.${seed}.RaPIDfiltered.map"
        genetic_map="${ibd_dir}/maps/chr${chr}.${pop}.${seed}.RaPIDfinal.map"

        ibd_out="${ibd_dir}/${pop}_chr${chr}_seed${seed}_out"
        ibd_bed="${ibd_dir}/${pop}_chr${chr}_seed${seed}_ibd.bed.gz"
        if [[ -s "$ibd_bed" ]]; then
                log_msg "${ibd_bed} already found for ${pop}; continue"
                continue
        fi

        # filtering gmap to this rapid run
        log_msg "filting gmap to this rapid run"
        log_msg "wc -l ${int_genetic_map} (before filter): $(wc -l < "$int_genetic_map")"
        python ~/software/RaPID/filter_mapping_file.py "$int_genetic_map" "$filt_genetic_map" 
        log_msg "wc -l ${filt_genetic_map} (after filter): $(wc -l < "$filt_genetic_map")"

        # interpolating gmap to this rapid run
        log_msg "interpolating gmap to this rapid run"
        n_vcf=$(bcftools view -H "$vcf" | wc -l)
        log_msg "wc -l ${vcf} (variants): ${n_vcf}"
        python ~/software/RaPID/interpolate_loci.py "$filt_genetic_map" "$vcf" "$genetic_map"
        log_msg "wc -l ${genetic_map} (final): $(wc -l < "$genetic_map")"

        # make sure genetic map lengths and vcf length are the same
        n_map=$(wc -l < "$genetic_map")
        if [[ "$n_vcf" -ne "$n_map" ]]; then
                log_msg "ERROR: VCF variant count (${n_vcf}) != genetic_map lines (${n_map}) for pop=${pop} seed=${seed}"
                exit 1
        fi

        # call ibd with rapid
        log_msg "Calling IBD of $pop for seed${seed} with rapid using...
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
                log_msg "no output for $pop seed$seed"
                exit 1
        fi

        # early exit if no ibd segments were called
        if [[ $(zcat "${ibd_out}/results.max.gz" | wc -l) -eq 0 ]]; then
                log_msg "no ibd segments detected for $pop on seed$seed; continue"
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

log_msg "called ibd for tmap=${tmap}, imap=${imap}, seed=${seed}"
