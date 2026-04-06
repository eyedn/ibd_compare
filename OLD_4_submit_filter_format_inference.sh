#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           submit_filter_format_inference.sh
###############################################################################

set -euo pipefail
source /home1/karatas/ibd_compare/other_scripts/const.sh
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh

# submit formatting for inference outputs
# ALL software outputs are formatted in this script
# NOTE: paths to software-param-specific dirs specifcied in this script
log_msg "Detected filter and format specifications
MIN_BP_LEN=${MIN_BP_LEN}
MIN_CM_LEN=${MIN_CM_LEN}
COMMON_FF_DIR=${COMMON_FF_DIR}
CHR=${CHR}"
TOTAL_TASKS=$(( NUM_REPS * ${#TMAPS[@]} * ${#IMAPS[@]} ))
sbatch \
    --job-name=infFF \
    --array="1-${TOTAL_TASKS}%${MAX_JOBS}" \
    job_scripts/filter_format_inference.sh \
        "$OUTDIR" \
        "$VCF_DIR" \
        "$COMMON_FF_DIR" \
        "$REF_DIR" \
        "$CHR" \
        "$MIN_BP_LEN" \
        "$MIN_CM_LEN" \
        "$HOT_BIN_SIZE" \
        "$HOT_BIN_SD" \
        "$NUM_REPS" \
        "${TMAPS[@]}" \
        -- \
        "${IMAPS[@]}" \
        -- \
        "${POPS[@]}"
log_msg "submitted formatting jobs for inferred ibd; TOTAL_TASKS = ${TOTAL_TASKS}"