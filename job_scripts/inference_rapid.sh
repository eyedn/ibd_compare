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
MIN_LEN="$1"
WINDOW_SIZE="$2"
NUM_RUNS="$3"
NUM_SUCCESSES="$4"
OUTDIR="$5"
VCF_DIR="$6"
NUM_REPS="$7"
CHR="$8"
shift 8

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

# define intermediate (chr, bp, cM) that will be filtered and interpolated
if [[ "$imap" == "HapMapII_GRCh38" ]]; then
        INT_GENETIC_MAP="/home1/karatas/references/RaPID_format_maps/HapMapII_GRCh38/plink.chr${CHR}.GRCh38.UCSC.reformated.map"
elif [[ "$imap" == "DeCodeSexAveraged_GRCh38" ]]; then
        INT_GENETIC_MAP="/home1/karatas/references/RaPID_format_maps/deCODE/chr${CHR}.reformated.map"
elif [[ "$imap" == "Flat1Mb1cM" ]]; then
        INT_GENETIC_MAP="/home1/karatas/references/RaPID_format_maps/1Mb1cM/chr${CHR}.reformated.map"
else
        echo "ERROR: unknown imap=${imap}" >&2
        exit 1
fi

# for each pop, infer ibd
ibd_dir="${OUTDIR}/rapid_sim_tmap${tmap}_imap${imap}_d${MIN_LEN}_w${WINDOW_SIZE}_r${NUM_RUNS}_s${NUM_SUCCESSES}"
mkdir -p "$ibd_dir"
mkdir -p "$ibd_dir/maps"
log_msg "ibd_dir defined as ${ibd_dir}"
for pop in "${POPS[@]}"; do
        vcf="${VCF_DIR}/${tmap}_${seed}_${pop}.vcf.gz"
        [[ -s "$vcf" ]] || { log_msg "vcf not found: $vcf"; exit 1; }

        # define filtered map and final genetic map paths after rapid filtering
        # and interpolating for this run
        filt_genetic_map="${ibd_dir}/maps/chr${CHR}.${pop}.${seed}.RaPIDfiltered.map"
        genetic_map="${ibd_dir}/maps/chr${CHR}.${pop}.${seed}.RaPIDfinal.map"

        ibd_out="${ibd_dir}/${pop}_chr${CHR}_seed${seed}_out"
        ibd_bed="${ibd_dir}/${pop}_chr${CHR}_seed${seed}_ibd.bed.gz"
        if [[ -s "$ibd_bed" ]]; then
                log_msg "${ibd_bed} already found for ${pop}; continue"
                continue
        fi

        # filtering gmap to this rapid run
        log_msg "filting gmap to this rapid run"
        log_msg "wc -l ${INT_GENETIC_MAP} (before filter): $(wc -l < "$INT_GENETIC_MAP")"
        python ~/software/RaPID/filter_mapping_file.py "$INT_GENETIC_MAP" "$filt_genetic_map" 
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
        min_length_in_cM = ${MIN_LEN}
        window_size = ${WINDOW_SIZE}
        num_runs = ${NUM_RUNS}
        num_successes = ${NUM_SUCCESSES}
        input_file_vcf_compressed = ${vcf}
        genetic_mapping_file = ${genetic_map}
        output_folder = $ibd_out"

        rm -rf "$ibd_out"
        "${HOME}/software/RaPID/RaPID_v.1.7" \
                -i "$vcf" \
                -g "$genetic_map" \
                -d "$MIN_LEN" \
                -w "$WINDOW_SIZE" \
                -r "$NUM_RUNS" \
                -s "$NUM_SUCCESSES" \
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
        TMPDIR="${SCRATCH:-/tmp}/bedtmp_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID:-$$}"
        mkdir -p "$TMPDIR"
        tmp_bed="${ibd_bed}.tmp"
        zcat "${ibd_out}/results.max.gz" | \
        awk -F'\t' 'BEGIN{OFS="\t"}   
                $6 ~ /^-?[0-9]+$/ && $7 ~ /^-?[0-9]+$/ && $6 >= 0 && $7 >= $6 {
                print $1,$6,$7,$2,$4,$3,$5,$8,$9,$10
                }' | \
        sort -k1,1 -k2,2n -k3,3n -T "$TMPDIR" -S 2G | \
        gzip > "$tmp_bed"

        mv -f "$tmp_bed" "$ibd_bed"
        rm -rf "$TMPDIR"
done

log_msg "called ibd for tmap=${tmap}, imap=${imap}, seed=${seed}"
