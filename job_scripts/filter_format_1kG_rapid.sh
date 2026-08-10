#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           filter_format_1kG_rapid.sh
###############################################################################

#SBATCH --time=1-00:00:00
#SBATCH --partition=qcb
#SBATCH --account=jazlynmo_738
#SBATCH --nodes=1
#SBATCH --output=/home1/karatas/logs/1kGrapFF/%A_%a.%x.out
#SBATCH --error=/home1/karatas/logs/1kGrapFF/%A_%a.%x.err
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
rapid_map_dir="$3"
common_ff_dir="$4"
ref_dir="$5"
min_bp_len="$6"
min_cm_len="$7"
hot_bin_size="$8"
hot_bin_sd="$9"
rap_min_len="${10}"
rap_window_size="${11}"
rap_num_runs="${12}"
rap_num_successes="${13}"
shift 13

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

rapid_dir_base="rapid_1kG_imap${imap}_d${rap_min_len}_w${rap_window_size}_r${rap_num_runs}_s${rap_num_successes}"
rapid_dir="${outdir}/${rapid_dir_base}"
rapid_ff_dir="${common_ff_dir}/${rapid_dir_base}"

mkdir -p "$common_ff_dir" "$rapid_ff_dir"

log_msg "filter_format_1kG_rapid task
task_id=${tid}
chr=${chr}
imap=${imap}
imap_label=${imap_label}
rapid_dir=${rapid_dir}
rapid_ff_dir=${rapid_ff_dir}
pops=${pops[*]}
ibdends_xmx=${ibdends_xmx}
ibdends_threads=${ibdends_threads}"

write_empty_standard_bed() {
    printf '%s\n' $'chr\tbp.str\tbp.end\tbp.len\tcM.str\tcM.end\tcM.len\tsam1\thap1\tsam2\thap2' > "$1"
}

has_gzip_records() {
    gzip -cd "$1" | awk 'NR { found = 1; exit } END { exit !found }'
}

for pop in "${pops[@]}"; do
    vcf="${vcf_dir}/chr${chr}/${pop}.vcf.gz"
    rapid_gmap="${rapid_map_dir}/pop${pop}_chr${chr}_imap${imap}.map"
    rapid_nonbed_dir="${rapid_dir}/${pop}_chr${chr}_out"
    rapid_nonbed="${rapid_nonbed_dir}/results.max.gz"
    rapid_nonbed_hapibdformat="${rapid_nonbed_dir}/results.max.hapibdformat.gz"
    rapid_ibd="${rapid_dir}/${pop}_chr${chr}_ibd.bed.gz"

    rapid_ibdends_prefix="${rapid_dir}/${pop}_chr${chr}.ibdends"
    rapid_ibdends_ibd="${rapid_ibdends_prefix}_ibd.ibdends.bed.gz"
    rapid_hot_hapibdformat="${rapid_nonbed_dir}/results.max.hot.hapibdformat.gz"
    rapid_hot_ibdends_prefix="${rapid_dir}/${pop}_chr${chr}.hot.ibdends"
    rapid_hot_ibdends_ibd="${rapid_hot_ibdends_prefix}_ibd.ibdends.bed.gz"

    rapid_prefix="${rapid_ff_dir}/rapid_1kG_${imap_label}_${pop}_chr${chr}"
    rapid_out="${rapid_prefix}.bed"
    rapid_hot_out="${rapid_prefix}.hot.bed"
    rapid_ibdends_out="${rapid_prefix}.ibdends.bed"
    rapid_hotibdends_out="${rapid_prefix}.hotibdends.bed"
    rapid_ibdendshot_out="${rapid_prefix}.ibdendshot.bed"

    [[ -s "$vcf" ]] || { echo "ERROR: missing vcf: $vcf" >&2; exit 1; }
    [[ -s "$rapid_nonbed" ]] || {
        echo "ERROR: missing rapid results: $rapid_nonbed" >&2
        exit 1
    }
    [[ -s "$rapid_ibd" ]] || {
        echo "ERROR: missing rapid bed.gz: $rapid_ibd" >&2
        exit 1
    }
    [[ -s "$rapid_gmap" ]] || {
        echo "ERROR: missing rapid gmap: $rapid_gmap" >&2
        exit 1
    }

    module purge
    ml openjdk/25.0.2

    log_msg "rapid: generating ${rapid_nonbed_hapibdformat}
pop=${pop}
chr=${chr}
imap=${imap_label}"
    zcat "$rapid_nonbed" | \
        awk 'BEGIN{OFS="\t"} ($6 >= 0 && $7 >= 0 && $6 <= $7) { print $2, $4+1, $3, $5+1, $1, $6, $7, $8 }' | \
        gzip -c > "$rapid_nonbed_hapibdformat"

    if ! has_gzip_records "$rapid_nonbed_hapibdformat"; then
        write_empty_standard_bed "$rapid_out"
        write_empty_standard_bed "$rapid_hot_out"
        write_empty_standard_bed "$rapid_ibdends_out"
        write_empty_standard_bed "$rapid_hotibdends_out"
        write_empty_standard_bed "$rapid_ibdendshot_out"
        gzip -f "$rapid_out" "$rapid_hot_out" "$rapid_ibdends_out" \
            "$rapid_hotibdends_out" "$rapid_ibdendshot_out"
        log_msg "no rapid segments; wrote empty standard outputs for ${pop}"
        continue
    fi

    log_msg "rapid: filtering segments using ibd-ends
pop=${pop}
chr=${chr}
imap=${imap_label}"
    java -Xmx"${ibdends_xmx}" -jar ~/software/ibd-ends/ibd-ends.jar \
        gt="$vcf" \
        ibd="$rapid_nonbed_hapibdformat" \
        map="$hapibd_gmap" \
        min-maf="${ibdends_min_maf}" \
        nthreads="${ibdends_threads}" \
        out="$rapid_ibdends_prefix"

    log_msg "rapid: generating ${rapid_ibdends_ibd}"
    zcat "${rapid_ibdends_prefix}.ibd.gz" | \
        awk 'BEGIN{OFS="\t"} NR>1 {print $5,$10,$11,$1,$2,$3,$4,$12}' | \
        gzip -c > "$rapid_ibdends_ibd"

    module purge
    ml rstats/4.5.1

    log_msg "rapid: format+filter raw calls
pop=${pop}
chr=${chr}
imap=${imap_label}"
    Rscript "/home1/karatas/ibd_compare/job_scripts/filter_and_format.R" \
        "$min_bp_len" "$min_cm_len" "$hot_bin_size" "$hot_bin_sd" \
        "$rapid_ibd" "$rapid_gmap" \
        "$rapid_out" "$rapid_hot_out" \
        "FALSE" "TRUE"

    log_msg "rapid: converting hot-filtered calls to hapibd format
pop=${pop}
chr=${chr}
imap=${imap_label}"
    Rscript "/home1/karatas/ibd_compare/job_scripts/convert_std_bed_to_hapibd.R" \
        "$rapid_hot_out" "$rapid_hot_hapibdformat"

    log_msg "rapid: format+filter ibd-ends adjusted calls
pop=${pop}
chr=${chr}
imap=${imap_label}"
    Rscript "/home1/karatas/ibd_compare/job_scripts/filter_and_format.R" \
        "$min_bp_len" "$min_cm_len" "$hot_bin_size" "$hot_bin_sd" \
        "$rapid_ibdends_ibd" "$hapibd_gmap" \
        "$rapid_ibdends_out" "$rapid_ibdendshot_out" \
        "TRUE" "TRUE"

    if ! has_gzip_records "$rapid_hot_hapibdformat"; then
        write_empty_standard_bed "$rapid_hotibdends_out"
        gzip -f "$rapid_out" "$rapid_hot_out" "$rapid_ibdends_out" \
            "$rapid_hotibdends_out" "$rapid_ibdendshot_out"
        log_msg "no hot-filtered rapid segments; wrote empty hot ibd-ends output"
        continue
    fi

    module purge
    ml openjdk/25.0.2

    log_msg "rapid: filtering hot-filtered segments using ibd-ends
pop=${pop}
chr=${chr}
imap=${imap_label}"
    java -Xmx"${ibdends_xmx}" -jar ~/software/ibd-ends/ibd-ends.jar \
        gt="$vcf" \
        ibd="$rapid_hot_hapibdformat" \
        map="$hapibd_gmap" \
        min-maf="${ibdends_min_maf}" \
        nthreads="${ibdends_threads}" \
        out="$rapid_hot_ibdends_prefix"

    log_msg "rapid: generating ${rapid_hot_ibdends_ibd}"
    zcat "${rapid_hot_ibdends_prefix}.ibd.gz" | \
        awk 'BEGIN{OFS="\t"} NR>1 {print $5,$10,$11,$1,$2,$3,$4,$12}' | \
        gzip -c > "$rapid_hot_ibdends_ibd"

    module purge
    ml rstats/4.5.1

    log_msg "rapid: formatting hot-filtered ibd-ends adjusted calls
pop=${pop}
chr=${chr}
imap=${imap_label}"
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
done

log_msg "1kG rapid filter+format done for chr=${chr}, imap=${imap}"
