#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           inference_hapibd.sh
###############################################################################

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --time=1-00:00:00
#SBATCH --partition=qcb
#SBATCH --account=jazlynmo_738
#SBATCH --nodes=1
#SBATCH --output=/home1/karatas/logs/hapibd.%A_%a.out
#SBATCH --error=/home1/karatas/logs/hapibd.%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=karatas@usc.edu

set -euo pipefail
module purge
ml legacy/CentOS7 jdk/17.0.5
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh

: "${SLURM_ARRAY_TASK_ID:?ERROR: run as a Slurm array}"

# input
MIN_SEED="$1"
MAX_GAP="$2"
MIN_EXTEND="$3"
MIN_OUTPUT="$4"
MIN_MARKERS="$5"
MIN_MAC="$6"
OUTDIR="$7"
VCF_DIR="$8"
NUM_REPS="$9"
CHR="${10}"
shift 10

# read TMAPS until "--", then read IMAPS until "--", then read POPS
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

# decode task id to which tmap + imap + seed
tid="${SLURM_ARRAY_TASK_ID}"
tid0=$(( tid - 1 ))

seed_id=$(( tid0 % NUM_REPS + 1 ))
tmp_id=$(( tid0 / NUM_REPS ))

imap_idx=$(( tmp_id % ${#IMAPS[@]} ))
tmap_idx=$(( tmp_id / ${#IMAPS[@]} ))

if (( tmap_idx < 0 || tmap_idx >= ${#TMAPS[@]} )); then
        echo "ERROR: task ${tid} produced tmap_idx=${tmap_idx}, but n_tmaps=${#TMAPS[@]}" >&2
        exit 1
fi
if (( imap_idx < 0 || imap_idx >= ${#IMAPS[@]} )); then
        echo "ERROR: task ${tid} produced imap_idx=${imap_idx}, but n_imaps=${#IMAPS[@]}" >&2
        exit 1
fi

tmap="${TMAPS[$tmap_idx]}"
imap="${IMAPS[$imap_idx]}"
seed="${seed_id}"

# define genetic map
if [[ "$imap" == "HapMapII_GRCh38" ]]; then
        GENETIC_MAP="/home1/karatas/references/HapMapII_GRCh38/plink.chr${CHR}.GRCh38.UCSC.format.map"
elif [[ "$imap" == "DeCodeSexAveraged_GRCh38" ]]; then
        GENETIC_MAP="/home1/karatas/references/deCODE/addendum/sex-averaged_gmaps/sex_av_chr${CHR}_hg38.gmap"
elif [[ "$imap" == "Flat1Mb1cM" ]]; then
        GENETIC_MAP="/home1/karatas/references/GRCh38.p14/recomb_maps_1Mb1cM/chr${CHR}.map"
else
        echo "ERROR: unknown imap=${imap}" >&2
        exit 1
fi

# for each pop, infer ibd
ibd_dir="${OUTDIR}/hapibd_sim_tmap${tmap}_imap${imap}_seed${MIN_SEED}_gap${MAX_GAP}_extend${MIN_EXTEND}_out${MIN_OUTPUT}_mark${MIN_MARKERS}_mac${MIN_MAC}"
mkdir -p "$ibd_dir"
log_msg "ibd_dir defined as ${ibd_dir}"

for pop in "${POPS[@]}"; do
        vcf="${VCF_DIR}/${tmap}_${seed}_${pop}.vcf.gz"
        [[ -s "$vcf" ]] || { log_msg "vcf missing/empty: $vcf"; exit 1; }

        prefix="${ibd_dir}/${pop}_chr${CHR}_seed${seed}"
        ibd_gz="${prefix}.ibd.gz"
        ibd_bed="${prefix}_ibd.bed.gz"

        if [[ -s "$ibd_bed" ]]; then
                log_msg "${ibd_bed} already found for ${pop}; continue"
                continue
        fi

        # call ibd with hapibd
        log_msg "Calling IBD of ${pop} for seed${seed} with hapibd using...
        min-seed=${MIN_SEED}
        max-gap=${MAX_GAP}
        min-extend=${MIN_EXTEND}
        min-output=${MIN_OUTPUT}
        min-markers=${MIN_MARKERS}
        min-mac=${MIN_MAC}
        gt=${vcf}
        map=${GENETIC_MAP}
        out=${prefix}"
        
        # hap-ibd writes directly to out=<prefix>
        rm -f "${prefix}.log" "${prefix}.hbd.gz" "${prefix}.ibd.gz" "${prefix}.ibd.gz.tbi" 2>/dev/null || true
        java -Xmx14G -jar ~/software/hap-ibd/hap-ibd.jar \
                min-seed="$MIN_SEED" \
                max-gap="$MAX_GAP" \
                min-extend="$MIN_EXTEND" \
                min-output="$MIN_OUTPUT" \
                min-markers="$MIN_MARKERS" \
                min-mac="$MIN_MAC" \
                nthreads="${SLURM_CPUS_PER_TASK:-1}" \
                gt="$vcf" \
                map="$GENETIC_MAP" \
                out="$prefix"

        # exit nonzero if hap-ibd did not produce an output
        if [[ ! -s "$ibd_gz" ]]; then
                log_msg "no hap-ibd output for ${pop} seed${seed}"
                exit 1
        fi

        # continue if no ibd segments were called (header-only or empty)
        if [[ $(zcat "$ibd_gz" | wc -l) -eq 0 ]]; then
                log_msg "no ibd segments detected for ${pop} on seed${seed}; continue"
                continue
        fi

        # convert ibd calls to bed
        TMPDIR="${SCRATCH:-/tmp}/bedtmp_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID:-$$}"
        mkdir -p "$TMPDIR"
        tmp_bed="${ibd_bed}.tmp"

        zcat "$ibd_gz" | \
        awk -F'\t' 'BEGIN{OFS="\t"} {print $5,$6,$7,$1,$2,$3,$4,$8}' | \
        sort -k1,1 -k2,2n -k3,3n -T "$TMPDIR" -S 2G | \
        gzip > "$tmp_bed"

        mv -f "$tmp_bed" "$ibd_bed"
        rm -rf "$TMPDIR"
done

log_msg "called ibd for tmap=${tmap}, imap=${imap}, seed=${seed}"
