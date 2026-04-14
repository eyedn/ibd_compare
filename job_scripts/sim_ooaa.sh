#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           sim_ooaa.sh
###############################################################################

#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --partition=qcb
#SBATCH --account=jazlynmo_738
#SBATCH --nodes=1
#SBATCH --output=/home1/karatas/logs/sim.%A_%a.%x.out
#SBATCH --error=/home1/karatas/logs/sim.%A_%a.%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=karatas@usc.edu

set -euo pipefail
module purge
ml gcc/13.3.0 htslib/1.19.1 conda
source /apps/conda/miniforge3/25.3.0/etc/profile.d/conda.sh
conda activate stdpopsim
export PATH="/home1/karatas/.conda/envs/stdpopsim/bin:$PATH"
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh

# constants
: "${SLURM_ARRAY_TASK_ID:?ERROR: run as a Slurm array}"

# inputs
tree_dir="$1"
sample_size="$2"
num_reps="$3"
msprime_model="$4"
chr="$5"
generation_time="$6"
mutation_rate="$7"
t_af_years="$8"
t_ooa_years="$9"
t_eu0_years="${10}"
t_eg_years="${11}"
r_eu0="${12}"
r_eu="${13}"
r_af="${14}"
n_a="${15}"
n_af1="${16}"
n_b="${17}"
n_eu0="${18}"
m_af_b="${19}"
m_af_eu="${20}"
admixture_time="${21}"
afr_prop="${22}"
eur_prop="${23}"
admix_initial_size="${24}"
admix_growth_rate_1="${25}"
admix_growth_rate_2="${26}"
admix_growth_rate_2_change_time="${27}"
census_time_offset="${28}"
shift 28
# decode task id to which map + replicate
maps=( "$@" )
if [[ ${#maps[@]} -eq 0 ]]; then
    echo "ERROR: no GMAPs provided (after shifting args)" >&2
    exit 1
fi
tid="${SLURM_ARRAY_TASK_ID}"
map_idx=$(( (tid - 1) / num_reps ))      # 0..n_maps-1
rep_id=$(( (tid - 1) % num_reps + 1 ))   # 1..num_reps
if (( map_idx < 0 || map_idx >= ${#maps[@]} )); then
    echo "ERROR: task ${tid} produced map_idx=${map_idx}, but n_maps=${#maps[@]}" >&2
    exit 1
fi
gmap="${maps[$map_idx]}"
seed="${rep_id}"

# check if this task was already completed; if so, skip
out_ts="${tree_dir}/${gmap}_${seed}_all.ts"
if [[ -s "${out_ts}.tsz" ]]; then
    log_msg "outputs exist (${out_ts}.tsz); skipping"
    exit 0
fi

# define sampling specification
sample_spec=(
    "AFR:${sample_size}"
    "EUR:${sample_size}"
    "ADMIX:${sample_size}"
)

mkdir -p "${tree_dir}"

log_msg "stdpopsim sim (combined array)
task_id=${tid}
gmap=${gmap} (map_idx=${map_idx}/${#maps[@]})
seed=${seed}
tree_dir=${tree_dir}
sample_size=${sample_size}
model=${msprime_model}
admixture_time=${admixture_time}
afr_prop=${afr_prop}
eur_prop=${eur_prop}
generation_time=${generation_time}
mutation_rate=${mutation_rate}
t_af_years=${t_af_years}
t_ooa_years=${t_ooa_years}
t_eu0_years=${t_eu0_years}
t_eg_years=${t_eg_years}
r_eu0=${r_eu0}
r_eu=${r_eu}
r_af=${r_af}
n_a=${n_a}
n_af1=${n_af1}
n_b=${n_b}
n_eu0=${n_eu0}
m_af_b=${m_af_b}
m_af_eu=${m_af_eu}
admix_initial_size=${admix_initial_size}
admix_growth_rate_1=${admix_growth_rate_1}
admix_growth_rate_2=${admix_growth_rate_2}
admix_growth_rate_2_change_time=${admix_growth_rate_2_change_time}
sample_spec=${sample_spec[*]}"

python /home1/karatas/ibd_compare/job_scripts/sim_ooaa.py \
    --output "$out_ts" \
    --sample-size "$sample_size" \
    --seed "$seed" \
    --msprime-model "$msprime_model" \
    --chromosome "$chr" \
    --genetic-map "$gmap" \
    --generation-time "$generation_time" \
    --mutation-rate "$mutation_rate" \
    --t-af-years "$t_af_years" \
    --t-ooa-years "$t_ooa_years" \
    --t-eu0-years "$t_eu0_years" \
    --t-eg-years "$t_eg_years" \
    --r-eu0 "$r_eu0" \
    --r-eu "$r_eu" \
    --r-af "$r_af" \
    --n-a "$n_a" \
    --n-af1 "$n_af1" \
    --n-b "$n_b" \
    --n-eu0 "$n_eu0" \
    --m-af-b "$m_af_b" \
    --m-af-eu "$m_af_eu" \
    --admixture-time "$admixture_time" \
    --afr-proportion "$afr_prop" \
    --eur-proportion "$eur_prop" \
    --admix-initial-size "$admix_initial_size" \
    --admix-growth-rate-1 "$admix_growth_rate_1" \
    --admix-growth-rate-2 "$admix_growth_rate_2" \
    --admix-growth-rate-2-change-time "$admix_growth_rate_2_change_time" \
    --census-time-offset "$census_time_offset"

# compress tree sequence
log_msg "compress tree sequence"
tszip "$out_ts"
rm -f "$out_ts"

log_msg "sim done"
