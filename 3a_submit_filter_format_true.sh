#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           submit_filter_format_true.sh
###############################################################################

set -euo pipefail
source /home1/karatas/ibd_compare/other_scripts/const.sh
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh

# submit formatting for true outputs
log_msg "Detected filter and format specifications
MIN_BP_LEN=${MIN_BP_LEN}
MIN_CM_LEN=${MIN_CM_LEN}
COMMON_FF_DIR=${COMMON_FF_DIR}
CHR=${CHR}"
TOTAL_TASKS=$(( NUM_REPS * ${#TMAPS[@]} ))
sbatch \
    --chdir=/home1/karatas/ibd_compare \
    --job-name=trueFF \
    --array="1-${TOTAL_TASKS}%${MAX_JOBS}" \
    job_scripts/filter_format_true.sh \
        "$TRUE_DIR" \
        "$COMMON_FF_DIR" \
        "$REF_DIR" \
        "$CHR" \
        "$MIN_BP_LEN" \
        "$MIN_CM_LEN" \
        "$NUM_REPS" \
        "${TMAPS[@]}" \
        -- \
        "${POPS[@]}"
log_msg "submitted formatting jobs for true ibd"
