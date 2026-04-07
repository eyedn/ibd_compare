#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           subset_ts.sh
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
ml gcc/13.3.0 htslib/1.19.1 bcftools/1.19 conda
source /apps/conda/miniforge3/25.3.0/etc/profile.d/conda.sh
conda activate stdpopsim
export PATH="/home1/karatas/.conda/envs/stdpopsim/bin:$PATH"
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh

: "${SLURM_ARRAY_TASK_ID:?ERROR: run as a Slurm array}"
THREADS="${SLURM_CPUS_PER_TASK:-1}"

# inputs
TREE_DIR="$1"
VCF_DIR="$2"
POP_INFO_DIR="$3"
CHR="$4"
SAMPLE_SIZE="$5"
NUM_REPS="$6"
shift 6

# read MAPS until "--", then read POPS after
MAPS=()
while [[ $# -gt 0 && "$1" != "--" ]]; do
    MAPS+=( "$1" )
    shift
done
if [[ $# -eq 0 ]]; then
    echo "ERROR: expected '--' separator between MAPS and POPS" >&2
    exit 1
fi
shift 1
POPS=( "$@" )

if [[ ${#MAPS[@]} -eq 0 ]]; then
    echo "ERROR: no MAPS provided" >&2
    exit 1
fi
if [[ ${#POPS[@]} -eq 0 ]]; then
    echo "ERROR: no POPS provided" >&2
    exit 1
fi

tid="${SLURM_ARRAY_TASK_ID}"
map_idx=$(( (tid - 1) / NUM_REPS ))
rep_id=$(( (tid - 1) % NUM_REPS + 1 ))

if (( map_idx < 0 || map_idx >= ${#MAPS[@]} )); then
    echo "ERROR: task ${tid} produced map_idx=${map_idx}, but n_maps=${#MAPS[@]}" >&2
    exit 1
fi

gmap="${MAPS[$map_idx]}"
seed="${rep_id}"

# check ts file was created
in_ts="${TREE_DIR}/${gmap}_${seed}_all.ts"
if [[ -s "${in_ts}.tsz" ]]; then
    log_msg "found compressed ts (${in_ts}.tsz); unzipping"
    tsunzip "${in_ts}.tsz"
elif [[ -s "${in_ts}" ]]; then
    log_msg "found uncompressed ts (${in_ts})"
else
    echo "ERROR: missing input tree sequence (${in_ts} or ${in_ts}.tsz)" >&2
    exit 1
fi

mkdir -p "${VCF_DIR}"
mkdir -p "${POP_INFO_DIR}"
log_msg "subset_ts task
TASK_ID=${tid}
GMAP=${gmap} (map_idx=${map_idx}/${#MAPS[@]})
SEED=${seed}
IN_TS=${in_ts}
POPS=${POPS[*]}"
all_vcf="${VCF_DIR}/${gmap}_${seed}_all.vcf.gz"
tmp_vcf="${VCF_DIR}/${gmap}_${seed}_all.vcf"

log_msg "generate vcf of all samples from tree (python export, 1-based)"
python - <<'PY' "${in_ts}" "${CHR}" "${tmp_vcf}"
import tskit, sys
ts = tskit.load(sys.argv[1])
sim_chr = sys.argv[2]
with open(sys.argv[3], "w") as f:
    ts.write_vcf(f, contig_id=f"chr{sim_chr}", position_transform=lambda xs: [x + 1 for x in xs])
PY
bcftools norm --rm-dup all --threads "$THREADS" "$tmp_vcf" | \
    bcftools view --types snps --min-alleles 2 --max-alleles 2 --threads "$THREADS" \
    -Oz -o "${all_vcf}"
tabix -p vcf "${all_vcf}"
rm -f "${tmp_vcf}"

for idx in "${!POPS[@]}"; do
    pop="${POPS[$idx]}"
    start=$(( idx * SAMPLE_SIZE ))
    end=$(( (idx + 1) * SAMPLE_SIZE - 1 ))
    sample_ids="${POP_INFO_DIR}/${gmap}_${seed}_${pop}_samples.txt"
    out_vcf="${VCF_DIR}/${gmap}_${seed}_${pop}.vcf.gz"

    : > "${sample_ids}"
    for i in $(seq "${start}" "${end}"); do
        echo "tsk_${i}" >> "${sample_ids}"
    done

    log_msg "subsetting ${pop} samples tsk_${start}..tsk_${end} -> ${out_vcf}"
    bcftools view \
        -S "${sample_ids}" \
        --threads "$THREADS" \
        -Oz -o "${out_vcf}" \
        "${all_vcf}"
    tabix -p vcf "${out_vcf}"
done

log_msg "rezipping input tree sequence to ${in_ts}.tsz"
tszip "${in_ts}"
rm -f "${in_ts}"

log_msg "subset_ts done for ${gmap} seed=${seed}"
