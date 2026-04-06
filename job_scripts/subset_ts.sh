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
CHR="$3"
NUM_REPS="$4"
shift 4

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

# for each pop, create a pop-specific ts from ts from sim
mkdir -p "${VCF_DIR}"
log_msg "subset_ts task
TASK_ID=${tid}
GMAP=${gmap} (map_idx=${map_idx}/${#MAPS[@]})
SEED=${seed}
IN_TS=${in_ts}
POPS=${POPS[*]}"
for pop in "${POPS[@]}"; do
    # hard coded pop ids from demographic model used
    case "${pop}" in
    AFR) pop_id=0 ;;
    EUR) pop_id=1 ;;
    ASIA) pop_id=2 ;;
    ADMIX) pop_id=3 ;;
    *) echo "ERROR: unknown pop '${pop}' (expected AFR/EUR/ASIA/ADMIX)" >&2; exit 1 ;;
    esac

    out_ts="${TREE_DIR}/${gmap}_${seed}_${pop}.ts"
    log_msg "subsetting tree to ${pop} (pop_id=${pop_id}) -> ${out_ts}"

    python - <<'PY' "${in_ts}" "${out_ts}" "${pop_id}"
import sys, tskit
in_ts, out_ts, pop_id = sys.argv[1], sys.argv[2], int(sys.argv[3])
ts = tskit.load(in_ts)
keep = [u for u in ts.samples() if ts.node(u).population == pop_id]
if not keep:
    raise RuntimeError(f"No sample nodes found for pop_id={pop_id}")
ts.simplify(keep, keep_unary=True).dump(out_ts)
PY
    log_msg "wrote ${out_ts}"

    # store compressed vcfs for each pop  (VCF-compliant 1-based positions)
    log_msg "generate vcf of pop-specific tree (python export, 1-based)"
    tmp_vcf="${VCF_DIR}/${gmap}_${seed}_${pop}.vcf"

    python - <<'PY' "${out_ts}" "${CHR}" "${tmp_vcf}"
import tskit, sys
ts = tskit.load(sys.argv[1])
sim_chr = sys.argv[2]
with open(sys.argv[3], "w") as f:
    ts.write_vcf(f, contig_id=f"chr{sim_chr}", position_transform=lambda xs: [x + 1 for x in xs])
PY
    bcftools norm --rm-dup all --threads "$THREADS" "$tmp_vcf" | \
        bcftools view --types snps --min-alleles 2 --max-alleles 2 --threads "$THREADS" \
        -Oz -o "${VCF_DIR}/${gmap}_${seed}_${pop}.vcf.gz"
    tabix -p vcf "${VCF_DIR}/${gmap}_${seed}_${pop}.vcf.gz"
    rm -f "${tmp_vcf}"

    tszip "${out_ts}"
    rm -f "${out_ts}"
    log_msg "wrote ${out_ts}.tsz"
done

log_msg "rezipping input tree sequence to ${in_ts}.tsz"
tszip "${in_ts}"
rm -f "${in_ts}"

log_msg "subset_ts done for ${gmap} seed=${seed}"
