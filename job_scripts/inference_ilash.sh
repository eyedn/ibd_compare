#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           inference_ilash.sh
###############################################################################

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --time=1-00:00:00
#SBATCH --partition=qcb
#SBATCH --account=jazlynmo_738
#SBATCH --nodes=1
#SBATCH --output=/home1/karatas/logs/ilash.%A_%a.out
#SBATCH --error=/home1/karatas/logs/ilash.%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=karatas@usc.edu

set -euo pipefail
module purge
ml gcc/13.3.0
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh

: "${SLURM_ARRAY_TASK_ID:?ERROR: run as a Slurm array}"

# inputs
auto_slice="$1"
slice_size="$2"
step_size="$3"
perm_count="$4"
shingle_size="$5"
shingle_overlap="$6"
bucket_count="$7"
match_threshold="$8"
interest_threshold="$9"
min_length="${10}"
slice_length="${11}"
cm_overlap="${12}"
minhash_threshold="${13}"
outdir="${14}"
ped_dir="${15}"
ilash_map_dir="${16}"
num_reps="${17}"
chr="${18}"
shift 18

tmaps=()
while [[ $# -gt 0 && "$1" != "--" ]]; do
    tmaps+=( "$1" )
    shift
done
[[ $# -gt 0 && "$1" == "--" ]] || {
    echo "ERROR: expected '--' after tmaps" >&2
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

[[ ${#tmaps[@]} -gt 0 ]] || { echo "ERROR: no tmaps provided" >&2; exit 1; }
[[ ${#imaps[@]} -gt 0 ]] || { echo "ERROR: no imaps provided" >&2; exit 1; }
[[ ${#pops[@]} -gt 0 ]] || { echo "ERROR: no pops provided" >&2; exit 1; }

tid="${SLURM_ARRAY_TASK_ID}"
tid0=$(( tid - 1 ))

seed_id=$(( tid0 % num_reps + 1 ))
tmp_id=$(( tid0 / num_reps ))

imap_idx=$(( tmp_id % ${#imaps[@]} ))
tmap_idx=$(( tmp_id / ${#imaps[@]} ))

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

ibd_dir="${outdir}/ilash_sim_tmap${tmap}_imap${imap}_auto${auto_slice}_min${min_length}_slice${slice_length}_overlap${cm_overlap}"
mkdir -p "${ibd_dir}/configs"
log_msg "ibd_dir defined as ${ibd_dir}"

for pop in "${pops[@]}"; do
    ped="${ped_dir}/${tmap}_${seed}_${pop}.ped"
    ilash_map="${ilash_map_dir}/pop${pop}_chr${chr}_seed${seed}_tmap${tmap}_imap${imap}.map"
    [[ -s "${ped}" ]] || { log_msg "ped missing/empty: ${ped}"; exit 1; }
    [[ -s "${ilash_map}" ]] || {
        log_msg "ilash map missing/empty: ${ilash_map}"
        exit 1
    }

    prefix="${ibd_dir}/${pop}_chr${chr}_seed${seed}"
    config_file="${ibd_dir}/configs/${pop}_chr${chr}_seed${seed}.config"
    raw_match="${prefix}.match"
    ilash_ibd_gz="${prefix}.ibd.gz"
    ibd_bed="${prefix}_ibd.bed.gz"

    if [[ -s "${ibd_bed}" ]]; then
        log_msg "${ibd_bed} already found for ${pop}; continue"
        continue
    fi

    cat > "${config_file}" <<EOF
map ${ilash_map}
ped ${ped}
output ${raw_match}
slice_size ${slice_size}
step_size ${step_size}
perm_count ${perm_count}
shingle_size ${shingle_size}
shingle_overlap ${shingle_overlap}
bucket_count ${bucket_count}
max_thread ${SLURM_CPUS_PER_TASK:-1}
match_threshold ${match_threshold}
interest_threshold ${interest_threshold}
min_length ${min_length}
auto_slice ${auto_slice}
slice_length ${slice_length}
cm_overlap ${cm_overlap}
minhash_threshold ${minhash_threshold}
EOF

    log_msg "Calling IBD of ${pop} for seed${seed} with ilash using...
    ped=${ped}
    map=${ilash_map}
    output=${raw_match}
    config=${config_file}"

    rm -f "${raw_match}" "${ilash_ibd_gz}" "${ibd_bed}" 2>/dev/null || true
    /home1/karatas/software/iLASH/build/ilash "${config_file}"

    if [[ ! -s "${raw_match}" ]]; then
        log_msg "no ilash output for ${pop} seed${seed}"
        exit 1
    fi

    if [[ $(wc -l < "${raw_match}") -eq 0 ]]; then
        log_msg "no ilash segments detected for ${pop} on seed${seed}; continue"
        continue
    fi

    awk 'BEGIN{OFS="\t"} {
        sam1 = $2
        hap1 = sam1
        sub(/^.*_/, "", hap1)
        sub(/_[^_]+$/, "", sam1)
        sam2 = $4
        hap2 = sam2
        sub(/^.*_/, "", hap2)
        sub(/_[^_]+$/, "", sam2)
        print sam1, hap1 + 1, sam2, hap2 + 1, $5, $6, $7, $10
    }' "${raw_match}" | gzip -c > "${ilash_ibd_gz}"

    tmpdir="${SCRATCH:-/tmp}/bedtmp_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID:-$$}"
    mkdir -p "${tmpdir}"
    tmp_bed="${ibd_bed}.tmp"

    zcat "${ilash_ibd_gz}" | \
    awk -F'\t' 'BEGIN{OFS="\t"} {print $5,$6,$7,$1,$2,$3,$4,$8}' | \
    sort -k1,1 -k2,2n -k3,3n -T "${tmpdir}" -S 2G | \
    gzip > "${tmp_bed}"

    mv -f "${tmp_bed}" "${ibd_bed}"
    rm -rf "${tmpdir}"
done

log_msg "called ibd for tmap=${tmap}, imap=${imap}, seed=${seed}"
