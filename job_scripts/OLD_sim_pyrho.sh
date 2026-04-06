#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           sim_pyrho.sh
###############################################################################

#SBATCH --cpus-per-task=1
#SBATCH --mem=48G
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
ml gcc/13.3.0 htslib/1.19.1 conda
source /apps/conda/miniforge3/25.3.0/etc/profile.d/conda.sh
conda activate stdpopsim
export PATH="/home1/karatas/.conda/envs/stdpopsim/bin:$PATH"
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh

GMAP="$1"
OUTDIR="$2"
POP="$3"

mkdir -p "$OUTDIR"

: "${SLURM_ARRAY_TASK_ID:?ERROR: run as a Slurm array}"
SEED="${SLURM_ARRAY_TASK_ID}"
THREADS="${SLURM_CPUS_PER_TASK:-1}"
MSPRIME_MODEL="${MSPRIME_MODEL:-dtwf}"

# simulation
log_msg "pyrho sim
GMAP=${GMAP}
POP=${POP}
SEED=${SEED}
OUTDIR=${OUTDIR}"
stdpopsim -v -e msprime --msprime-model "${MSPRIME_MODEL}" \
  HomSap -s "${SEED}" -g "${GMAP}" -c chr1 \
  --left 0 --right 100000000 \
  -b "${OUTDIR}/${GMAP}_${SEED}_${POP}.bib" \
  -o "${OUTDIR}/${GMAP}_${SEED}_${POP}.ts" \
  -d OutOfAfricaArchaicAdmixture_5R19 \
  "${POP}:5000"

# store compressed vcf
log_msg "generate vcf for ${POP}"
tskit vcf --contig-id chr1 "${OUTDIR}/${GMAP}_${SEED}_${POP}.ts" \
  | bgzip -@ "${THREADS}" -c > "${OUTDIR}/${GMAP}_${SEED}_${POP}.vcf.gz"
tabix -p vcf "${OUTDIR}/${GMAP}_${SEED}_${POP}.vcf.gz"

# compress tree
log_msg "compress tree sequence"
tszip "${OUTDIR}/${GMAP}_${SEED}_${POP}.ts"
log_msg "Pyrho SIM ${POP} done"
