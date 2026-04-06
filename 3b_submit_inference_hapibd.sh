#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           submit_inference_hapibd.sh
###############################################################################

set -euo pipefail
source /home1/karatas/ibd_compare/other_scripts/const.sh
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh

# constants
MAX_JOBS=50
VCF_DIR="${OUTDIR}/vcfs"

# submit ibd inference
log_msg "Detected hapibd params
SIM_CHR=${CHR}
OUTDIR=${OUTDIR}
MIN_SEED=${HAP_MIN_SEED} 
MAX_GAP=${HAP_MAX_GAP}
MIN_EXTEND=${HAP_MIN_EXTEND}
MIN_OUTPUT=${HAP_MIN_OUTPUT}
MIN_MARKERS=${HAP_MIN_MARKERS}
MIN_MAC=${HAP_MIN_MAC}"
TOTAL_TASKS=$(( NUM_REPS * ${#TMAPS[@]} * ${#IMAPS[@]} ))
sbatch \
  --chdir=/home1/karatas/ibd_compare \
  --job-name=hapibd \
  --array="1-${TOTAL_TASKS}%${MAX_JOBS}" \
  --exclude=d11-05,d11-06,d11-07,b23-11,b23-13 \
  job_scripts/inference_hapibd.sh \
    "$HAP_MIN_SEED" \
    "$HAP_MAX_GAP" \
    "$HAP_MIN_EXTEND" \
    "$HAP_MIN_OUTPUT" \
    "$HAP_MIN_MARKERS" \
    "$HAP_MIN_MAC" \
    "$OUTDIR" \
    "$VCF_DIR" \
    "$NUM_REPS" \
    "$CHR" \
    "${TMAPS[@]}" \
    -- \
    "${IMAPS[@]}" \
    -- \
    "${POPS[@]}"
log_msg "submitted hapibd inference TOTAL_TASKS=${TOTAL_TASKS}"
