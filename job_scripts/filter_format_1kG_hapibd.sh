#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           filter_format_1kG_hapibd.sh
###############################################################################

#SBATCH --time=1-00:00:00
#SBATCH --partition=qcb
#SBATCH --account=jazlynmo_738
#SBATCH --nodes=1
#SBATCH --output=/home1/karatas/logs/1kGhapFF/%A_%a.%x.out
#SBATCH --error=/home1/karatas/logs/1kGhapFF/%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=karatas@usc.edu

set -euo pipefail
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh

: "${IBDENDS_XMX:?ERROR: IBDENDS_XMX must be set by submit wrapper}"
: "${IBDENDS_THREADS:?ERROR: IBDENDS_THREADS must be set by submit wrapper}"
: "${IBDENDS_MIN_MAF:?ERROR: IBDENDS_MIN_MAF must be set by submit wrapper}"
: "${SLURM_ARRAY_TASK_ID:?ERROR: run as a Slurm array}"

ibdends_xmx="${IBDENDS_XMX}"
ibdends_threads="${IBDENDS_THREADS}"
ibdends_min_maf="${IBDENDS_MIN_MAF}"

# inputs
outdir="$1"
vcf_dir="$2"
common_ff_dir="$3"
ref_dir="$4"
min_bp_len="$5"
min_cm_len="$6"
hot_bin_size="$7"
hot_bin_sd="$8"
hap_min_seed="$9"
hap_max_gap="${10}"
hap_min_extend="${11}"
hap_min_output="${12}"
hap_min_markers="${13}"
hap_min_mac="${14}"
shift 14

chrs=()
while [[ $# -gt 0 && "$1" != "--" ]]; do
    chrs+=( "$1" )
    shift
done
[[ $# -gt 0 && "$1" == "--" ]] || {
    echo "ERROR: expected '--' after chromosomes" >&2
    exit 1
}
shift

imaps=()
while [[ $# -gt 0 && "$1" != "--" ]]; do
    imaps+=( "$1" )
    shift
done
[[ $# -gt 0 && "$1" == "--" ]] || {
    echo "ERROR: expected '--' after imaps" >&2
    exit 1
}
shift

pops=( "$@" )

[[ ${#chrs[@]} -gt 0 ]] || { echo "ERROR: no chromosomes provided" >&2; exit 1; }
[[ ${#imaps[@]} -gt 0 ]] || { echo "ERROR: no imaps provided" >&2; exit 1; }
[[ ${#pops[@]} -gt 0 ]] || { echo "ERROR: no pops provided" >&2; exit 1; }

# decode task id -> (chr, imap)
tid="${SLURM_ARRAY_TASK_ID}"
tid0=$(( tid - 1 ))
imap_idx=$(( tid0 % ${#imaps[@]} ))
chr_idx=$(( tid0 / ${#imaps[@]} ))

if (( chr_idx < 0 || chr_idx >= ${#chrs[@]} )); then
    echo "ERROR: task ${tid} produced chr_idx=${chr_idx}" >&2
    exit 1
fi
if (( imap_idx < 0 || imap_idx >= ${#imaps[@]} )); then
    echo "ERROR: task ${tid} produced imap_idx=${imap_idx}" >&2
    exit 1
fi

chr="${chrs[$chr_idx]}"
imap="${imaps[$imap_idx]}"

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
    echo "ERROR: unknown imap=${imap}" >&2
    exit 1
fi

[[ -s "$hapibd_gmap" ]] || {
    echo "ERROR: missing hapibd_gmap: $hapibd_gmap" >&2
    exit 1
}

hapibd_dir_base="hapibd_1kG_imap${imap}_seed${hap_min_seed}_gap${hap_max_gap}_extend${hap_min_extend}_out${hap_min_output}_mark${hap_min_markers}_mac${hap_min_mac}"
hapibd_dir="${outdir}/${hapibd_dir_base}"
hapibd_ff_dir="${common_ff_dir}/${hapibd_dir_base}"

mkdir -p "$common_ff_dir" "$hapibd_ff_dir"

# Write the standard BED schema for paths with no callable segments.
write_empty_standard_bed() {
    printf '%s\n' \
        $'chr\tbp.str\tbp.end\tbp.len\tcM.str\tcM.end\tcM.len\tsam1\thap1\tsam2\thap2' \
        > "$1"
}

log_msg "filter_format_1kG_hapibd task
task_id=${tid}
chr=${chr}
imap=${imap}
imap_label=${imap_label}
hapibd_dir=${hapibd_dir}
hapibd_ff_dir=${hapibd_ff_dir}
pops=${pops[*]}
ibdends_xmx=${ibdends_xmx}
ibdends_threads=${ibdends_threads}"

for pop in "${pops[@]}"; do
    vcf="${vcf_dir}/chr${chr}/${pop}.vcf.gz"
    hapibd_nonbed="${hapibd_dir}/${pop}_chr${chr}.ibd.gz"
    hapibd_ibd="${hapibd_dir}/${pop}_chr${chr}_ibd.bed.gz"

    hapibd_ibdends_prefix="${hapibd_dir}/${pop}_chr${chr}.ibdends"
    hapibd_ibdends_ibd="${hapibd_ibdends_prefix}_ibd.ibdends.bed.gz"
    hapibd_hot_hapibdformat="${hapibd_dir}/${pop}_chr${chr}.hot.hapibdformat.gz"
    hapibd_hot_ibdends_prefix="${hapibd_dir}/${pop}_chr${chr}.hot.ibdends"
    hapibd_hot_ibdends_ibd="${hapibd_hot_ibdends_prefix}_ibd.ibdends.bed.gz"

    hapibd_prefix="${hapibd_ff_dir}/hapibd_1kG_${imap_label}_${pop}_chr${chr}"
    hapibd_out="${hapibd_prefix}.bed"
    hapibd_hot_out="${hapibd_prefix}.hot.bed"
    hapibd_ibdends_out="${hapibd_prefix}.ibdends.bed"
    hapibd_hotibdends_out="${hapibd_prefix}.hotibdends.bed"
    hapibd_ibdendshot_out="${hapibd_prefix}.ibdendshot.bed"

    [[ -s "$vcf" ]] || { echo "ERROR: missing vcf: $vcf" >&2; exit 1; }
    [[ -s "$hapibd_nonbed" ]] || {
        echo "ERROR: missing hapibd ibd.gz: $hapibd_nonbed" >&2
        exit 1
    }
    [[ -s "$hapibd_ibd" ]] || {
        echo "ERROR: missing hapibd bed.gz: $hapibd_ibd" >&2
        exit 1
    }

    # Format raw calls before ibd-ends so an empty raw result can short-circuit.
    module purge
    ml rstats/4.5.1
    log_msg "hapibd: format+filter raw calls
pop=${pop}
chr=${chr}
imap=${imap_label}"
    Rscript "/home1/karatas/ibd_compare/job_scripts/filter_and_format.R" \
        "$min_bp_len" "$min_cm_len" "$hot_bin_size" "$hot_bin_sd" \
        "$hapibd_ibd" "$hapibd_gmap" \
        "$hapibd_out" "$hapibd_hot_out" \
        "TRUE" "TRUE"

    if [[ $(tail -n +2 "$hapibd_out" | wc -l) -eq 0 ]]; then
        log_msg "hapibd: raw calls empty; skipping ibd-ends"
        write_empty_standard_bed "$hapibd_ibdends_out"
        write_empty_standard_bed "$hapibd_hotibdends_out"
        write_empty_standard_bed "$hapibd_ibdendshot_out"
        gzip -f "$hapibd_out" "$hapibd_hot_out" "$hapibd_ibdends_out" \
            "$hapibd_hotibdends_out" "$hapibd_ibdendshot_out"
        continue
    fi

    module purge
    ml openjdk/25.0.2

    log_msg "hapibd: filtering segments using ibd-ends
pop=${pop}
chr=${chr}
imap=${imap_label}"
    java -Xmx"${ibdends_xmx}" -jar ~/software/ibd-ends/ibd-ends.jar \
        gt="$vcf" \
        ibd="$hapibd_nonbed" \
        map="$hapibd_gmap" \
        min-maf="${ibdends_min_maf}" \
        nthreads="${ibdends_threads}" \
        out="$hapibd_ibdends_prefix"

    log_msg "hapibd: generating ${hapibd_ibdends_ibd}"
    zcat "${hapibd_ibdends_prefix}.ibd.gz" | \
        awk 'BEGIN{OFS="\t"} NR>1 {print $5,$10,$11,$1,$2,$3,$4,$12}' | \
        gzip -c > "$hapibd_ibdends_ibd"

    module purge
    ml rstats/4.5.1
    log_msg "hapibd: converting hot-filtered calls to hapibd format
pop=${pop}
chr=${chr}
imap=${imap_label}"
    Rscript "/home1/karatas/ibd_compare/job_scripts/convert_std_bed_to_hapibd.R" \
        "$hapibd_hot_out" "$hapibd_hot_hapibdformat"

    if [[ $(tail -n +2 "$hapibd_hot_out" | wc -l) -eq 0 ]]; then
        log_msg "hapibd: hot-filtered calls empty; skipping hot ibd-ends"
        write_empty_standard_bed "$hapibd_hotibdends_out"
        gzip -f "$hapibd_out" "$hapibd_hot_out" "$hapibd_ibdends_out" \
            "$hapibd_hotibdends_out" "$hapibd_ibdendshot_out"
        continue
    fi

    module purge
    ml openjdk/25.0.2

    log_msg "hapibd: filtering hot-filtered segments using ibd-ends
pop=${pop}
chr=${chr}
imap=${imap_label}"
    java -Xmx"${ibdends_xmx}" -jar ~/software/ibd-ends/ibd-ends.jar \
        gt="$vcf" \
        ibd="$hapibd_hot_hapibdformat" \
        map="$hapibd_gmap" \
        min-maf="${ibdends_min_maf}" \
        nthreads="${ibdends_threads}" \
        out="$hapibd_hot_ibdends_prefix"

    log_msg "hapibd: generating ${hapibd_hot_ibdends_ibd}"
    zcat "${hapibd_hot_ibdends_prefix}.ibd.gz" | \
        awk 'BEGIN{OFS="\t"} NR>1 {print $5,$10,$11,$1,$2,$3,$4,$12}' | \
        gzip -c > "$hapibd_hot_ibdends_ibd"

    module purge
    ml rstats/4.5.1

    log_msg "hapibd: format+filter ibd-ends adjusted calls
pop=${pop}
chr=${chr}
imap=${imap_label}"
    Rscript "/home1/karatas/ibd_compare/job_scripts/filter_and_format.R" \
        "$min_bp_len" "$min_cm_len" "$hot_bin_size" "$hot_bin_sd" \
        "$hapibd_ibdends_ibd" "$hapibd_gmap" \
        "$hapibd_ibdends_out" "$hapibd_ibdendshot_out" \
        "TRUE" "TRUE"

    log_msg "hapibd: formatting hot-filtered ibd-ends adjusted calls
pop=${pop}
chr=${chr}
imap=${imap_label}"
    Rscript "/home1/karatas/ibd_compare/job_scripts/filter_and_format.R" \
        "$min_bp_len" "$min_cm_len" "$hot_bin_size" "$hot_bin_sd" \
        "$hapibd_hot_ibdends_ibd" "$hapibd_gmap" \
        "$hapibd_hotibdends_out" "/tmp/hapibd_hotibdends_unused.bed" \
        "TRUE" "FALSE"

    gzip -f \
        "$hapibd_out" \
        "$hapibd_hot_out" \
        "$hapibd_ibdends_out" \
        "$hapibd_hotibdends_out" \
        "$hapibd_ibdendshot_out"
done

log_msg "1kG hapibd filter+format done for chr=${chr}, imap=${imap}"
