#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           prepare_pop_vcfs_peds.sh
###############################################################################

#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=1-00:00:00
#SBATCH --partition=qcb
#SBATCH --account=jazlynmo_738
#SBATCH --nodes=1
#SBATCH --output=/home1/karatas/logs/slurm.%A_%a.%x.out
#SBATCH --error=/home1/karatas/logs/slurm.%A_%a.%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=karatas@usc.edu

set -euo pipefail
module purge
ml gcc/13.3.0 htslib/1.19.1 bcftools/1.19 conda plink2/2.00a4.3
source /apps/conda/miniforge3/25.3.0/etc/profile.d/conda.sh
conda activate stdpopsim
export PATH="/home1/karatas/.conda/envs/stdpopsim/bin:$PATH"
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh

: "${SLURM_ARRAY_TASK_ID:?ERROR: run as a Slurm array}"
threads="${SLURM_CPUS_PER_TASK:-1}"

# inputs
tree_dir="$1"
vcf_dir="$2"
ped_dir="$3"
pop_info_dir="$4"
chr="$5"
sample_size="$6"
num_reps="$7"
shift 7

# read maps until "--", then read pops after
maps=()
while [[ $# -gt 0 && "$1" != "--" ]]; do
    maps+=( "$1" )
    shift
done
if [[ $# -eq 0 ]]; then
    echo "ERROR: expected '--' separator between maps and pops" >&2
    exit 1
fi
shift 1
pops=( "$@" )

if [[ ${#maps[@]} -eq 0 ]]; then
    echo "ERROR: no maps provided" >&2
    exit 1
fi
if [[ ${#pops[@]} -eq 0 ]]; then
    echo "ERROR: no pops provided" >&2
    exit 1
fi

tid="${SLURM_ARRAY_TASK_ID}"
map_idx=$(( (tid - 1) / num_reps ))
rep_id=$(( (tid - 1) % num_reps + 1 ))

if (( map_idx < 0 || map_idx >= ${#maps[@]} )); then
    echo "ERROR: task ${tid} produced map_idx=${map_idx}, but n_maps=${#maps[@]}" >&2
    exit 1
fi

gmap="${maps[$map_idx]}"
seed="${rep_id}"

need_tree=false
for idx in "${!pops[@]}"; do
    pop="${pops[$idx]}"
    out_vcf="${vcf_dir}/${gmap}_${seed}_${pop}.vcf.gz"
    if [[ ! -s "${out_vcf}" ]]; then
        need_tree=true
        break
    fi
done

mkdir -p "${vcf_dir}"
mkdir -p "${ped_dir}"
mkdir -p "${pop_info_dir}"

all_vcf="${vcf_dir}/${gmap}_${seed}_all.vcf.gz"
tmp_vcf="${vcf_dir}/${gmap}_${seed}_all.vcf"

temp_dir=""
temp_ts=""

if [[ "${need_tree}" == "true" ]]; then
    # stage a task-local copy of the tree sequence so shared inputs remain
    # read-only across downstream jobs
    in_ts="${tree_dir}/${gmap}_${seed}_all.ts"
    temp_root="${TMPDIR:-/tmp}"
    temp_dir="$(mktemp -d \
        "${temp_root}/prepare_pop_vcfs_peds.${gmap}_${seed}_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.XXXXXX")"
    temp_ts="${temp_dir}/${gmap}_${seed}_all.ts"
    cleanup() {
        rm -rf "${temp_dir}"
    }
    trap cleanup EXIT

    if [[ -s "${in_ts}.tsz" ]]; then
        temp_tsz="${temp_dir}/${gmap}_${seed}_all.ts.tsz"
        cp "${in_ts}.tsz" "${temp_tsz}"
        log_msg "staging compressed ts
canonical_input=${in_ts}.tsz
temp_input=${temp_tsz}
temp_output=${temp_ts}"
        tsunzip "${temp_tsz}"
    elif [[ -s "${in_ts}" ]]; then
        log_msg "staging uncompressed ts
canonical_input=${in_ts}
temp_output=${temp_ts}"
        cp "${in_ts}" "${temp_ts}"
    else
        echo "ERROR: missing input tree sequence (${in_ts} or ${in_ts}.tsz)" >&2
        exit 1
    fi
fi

log_msg "prepare_pop_vcfs_peds task
task_id=${tid}
gmap=${gmap} (map_idx=${map_idx}/${#maps[@]})
seed=${seed}
pops=${pops[*]}"

if [[ "${need_tree}" == "true" && ! -s "${all_vcf}" ]]; then
log_msg "generate vcf of all samples from tree (python export, 1-based)"
python - <<'PY' "${temp_ts}" "${chr}" "${tmp_vcf}"
import tskit, sys
ts = tskit.load(sys.argv[1])
sim_chr = sys.argv[2]
with open(sys.argv[3], "w") as f:
    ts.write_vcf(f, contig_id=f"chr{sim_chr}", position_transform=lambda xs: [x + 1 for x in xs])
PY
    bcftools norm --rm-dup all --threads "$threads" "$tmp_vcf" | \
        bcftools view --types snps --min-alleles 2 --max-alleles 2 \
        --threads "$threads" -Oz -o "${all_vcf}"
    tabix -p vcf "${all_vcf}"
    rm -f "${tmp_vcf}"
fi

for idx in "${!pops[@]}"; do
    pop="${pops[$idx]}"
    start=$(( idx * sample_size ))
    end=$(( (idx + 1) * sample_size - 1 ))
    sample_ids="${pop_info_dir}/${gmap}_${seed}_${pop}_samples.txt"
    out_vcf="${vcf_dir}/${gmap}_${seed}_${pop}.vcf.gz"
    ped_prefix="${ped_dir}/${gmap}_${seed}_${pop}"
    ped_file="${ped_prefix}.ped"
    ped_map="${ped_prefix}.map"

    if [[ ! -s "${sample_ids}" ]]; then
        : > "${sample_ids}"
        for i in $(seq "${start}" "${end}"); do
            echo "tsk_${i}" >> "${sample_ids}"
        done
    fi

    if [[ -s "${out_vcf}" ]]; then
        log_msg "vcf reused (${out_vcf})"
    else
        [[ -s "${all_vcf}" ]] || {
            echo "ERROR: missing all-sample VCF ${all_vcf}" >&2
            exit 1
        }
        log_msg "subsetting ${pop} samples tsk_${start}..tsk_${end} -> ${out_vcf}"
        bcftools view \
            -S "${sample_ids}" \
            --threads "$threads" \
            -Oz -o "${out_vcf}" \
            "${all_vcf}"
        tabix -p vcf "${out_vcf}"
        log_msg "vcf created (${out_vcf})"
    fi

    if [[ -s "${ped_file}" && -s "${ped_map}" ]]; then
        log_msg "ped reused (${ped_file})"
        continue
    fi

    log_msg "creating ped for ${gmap}_${seed}_${pop}"
    plink2 \
        --vcf "${out_vcf}" \
        --export ped \
        --out "${ped_prefix}"
    [[ -s "${ped_file}" && -s "${ped_map}" ]] || {
        echo "ERROR: expected PED outputs missing for ${ped_prefix}" >&2
        exit 1
    }
    log_msg "ped created (${ped_file})"
done

log_msg "prepare_pop_vcfs_peds done for ${gmap} seed=${seed}"
