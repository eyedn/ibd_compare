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
#SBATCH --output=/home1/karatas/logs/hapibd/%A_%a.%x.out
#SBATCH --error=/home1/karatas/logs/hapibd/%A_%a.%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=karatas@usc.edu

set -euo pipefail
module purge
ml openjdk/25.0.2
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh

: "${SLURM_ARRAY_TASK_ID:?ERROR: run as a Slurm array}"

# input
min_seed="$1"
max_gap="$2"
min_extend="$3"
min_output="$4"
min_markers="$5"
min_mac="$6"
outdir="$7"
vcf_dir="$8"
num_reps="$9"
chr="${10}"
contract_file="${11}"
id_contract="${12}"
shift 12

source /home1/karatas/ibd_compare/other_scripts/require_sim_id_contract.sh \
    "$contract_file" "$id_contract"

# read tmaps until "--", then read imaps until "--", then read pops
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
        echo "ERROR: task ${tid} produced tmap_idx=${tmap_idx}, but n_tmaps=${#tmaps[@]}" >&2
        exit 1
fi
if (( imap_idx < 0 || imap_idx >= ${#imaps[@]} )); then
        echo "ERROR: task ${tid} produced imap_idx=${imap_idx}, but n_imaps=${#imaps[@]}" >&2
        exit 1
fi

tmap="${tmaps[$tmap_idx]}"
imap="${imaps[$imap_idx]}"

# define genetic map
if [[ "$imap" == "HapMapII_GRCh38" ]]; then
        genetic_map="/home1/karatas/references/HapMapII_GRCh38/plink.chr${chr}.GRCh38.UCSC.format.map"
elif [[ "$imap" == "DeCodeSexAveraged_GRCh38" ]]; then
        genetic_map="/home1/karatas/references/deCODE/addendum/sex-averaged_gmaps/sex_av_chr${chr}_hg38.gmap"
elif [[ "$imap" == "Flat1Mb1cM" ]]; then
        genetic_map="/home1/karatas/references/GRCh38.p14/recomb_maps_1Mb1cM/chr${chr}.map"
else
        echo "ERROR: unknown imap=${imap}" >&2
        exit 1
fi

# for each pop, infer ibd
ibd_dir="${outdir}/hapibd_sim_tmap${tmap}_imap${imap}_minseed${min_seed}_gap${max_gap}_extend${min_extend}_out${min_output}_mark${min_markers}_mac${min_mac}"
mkdir -p "$ibd_dir"
log_msg "ibd_dir defined as ${ibd_dir}"

for pop in "${pops[@]}"; do
        vcf="${vcf_dir}/${tmap}_${rep}_chr${chr}_${pop}.vcf.gz"
        [[ -s "$vcf" ]] || { log_msg "vcf missing/empty: $vcf"; exit 1; }

        prefix="${ibd_dir}/${pop}_chr${chr}_rep${rep}"
        ibd_gz="${prefix}.ibd.gz"
        ibd_bed="${prefix}_ibd.bed.gz"

        if [[ -s "$ibd_bed" ]]; then
                log_msg "${ibd_bed} already found for ${pop}; continue"
                continue
        fi

        # call ibd with hapibd
        log_msg "Calling IBD of ${pop} for rep${rep} with hapibd using...
        min-seed=${min_seed}
        max-gap=${max_gap}
        min-extend=${min_extend}
        min-output=${min_output}
        min-markers=${min_markers}
        min-mac=${min_mac}
        gt=${vcf}
        map=${genetic_map}
        out=${prefix}"
        
        # hap-ibd writes directly to out=<prefix>
        rm -f "${prefix}.log" "${prefix}.hbd.gz" "${prefix}.ibd.gz" "${prefix}.ibd.gz.tbi" 2>/dev/null || true
        java -Xmx14G -jar ~/software/hap-ibd/hap-ibd.jar \
                min-seed="$min_seed" \
                max-gap="$max_gap" \
                min-extend="$min_extend" \
                min-output="$min_output" \
                min-markers="$min_markers" \
                min-mac="$min_mac" \
                nthreads="${SLURM_CPUS_PER_TASK:-1}" \
                gt="$vcf" \
                map="$genetic_map" \
                out="$prefix"

        # exit nonzero if hap-ibd did not produce an output
        if [[ ! -s "$ibd_gz" ]]; then
                log_msg "no hap-ibd output for ${pop} rep${rep}"
                exit 1
        fi

        # Preserve an empty BED artifact for downstream format/filter workers.
        if [[ $(zcat "$ibd_gz" | wc -l) -eq 0 ]]; then
                log_msg "no ibd segments detected for ${pop} on rep${rep}"
                gzip -c /dev/null > "$ibd_bed"
                continue
        fi

        # convert ibd calls to bed
        tmpdir="${SCRATCH:-/tmp}/bedtmp_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID:-$$}"
        mkdir -p "$tmpdir"
        tmp_bed="${ibd_bed}.tmp"

        zcat "$ibd_gz" | \
        awk -F'\t' 'BEGIN{OFS="\t"} {print $5,$6,$7,$1,$2,$3,$4,$8}' | \
        sort -k1,1 -k2,2n -k3,3n -T "$tmpdir" -S 2G | \
        gzip > "$tmp_bed"

        mv -f "$tmp_bed" "$ibd_bed"
        rm -rf "$tmpdir"
done

log_msg "called ibd for tmap=${tmap}, imap=${imap}, rep=${rep}"
