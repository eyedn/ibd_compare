#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           2b_submit_inference_hapibd.sh
###############################################################################


set -euo pipefail
source /home1/karatas/ibd_compare/other_scripts/const.sh
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh

# constants
max_jobs=50
vcf_dir="${OUTDIR}/vcfs"

# submit ibd inference
log_msg "Detected hapibd params
sim_chr=${CHR}
outdir=${OUTDIR}
min_seed=${HAP_MIN_SEED} 
max_gap=${HAP_MAX_GAP}
min_extend=${HAP_MIN_EXTEND}
min_output=${HAP_MIN_OUTPUT}
min_markers=${HAP_MIN_MARKERS}
min_mac=${HAP_MIN_MAC}"
total_tasks=$(( NUM_REPS * ${#TMAPS[@]} * ${#IMAPS[@]} ))
sbatch \
  --chdir=/home1/karatas/ibd_compare \
  --job-name=hapibd \
  --array="1-${total_tasks}%${max_jobs}" \
  --exclude=d11-05,d11-06,d11-07,b23-11,b23-13 \
  job_scripts/inference_hapibd.sh \
    "$HAP_MIN_SEED" \
    "$HAP_MAX_GAP" \
    "$HAP_MIN_EXTEND" \
    "$HAP_MIN_OUTPUT" \
    "$HAP_MIN_MARKERS" \
    "$HAP_MIN_MAC" \
    "$OUTDIR" \
    "$vcf_dir" \
    "$NUM_REPS" \
    "$CHR" \
    "$SIM_ID_CONTRACT_FILE" \
    "$SIM_ID_CONTRACT" \
    "${TMAPS[@]}" \
    -- \
    "${IMAPS[@]}" \
    -- \
    "${POPS[@]}"
log_msg "submitted hapibd inference total_tasks=${total_tasks}"
