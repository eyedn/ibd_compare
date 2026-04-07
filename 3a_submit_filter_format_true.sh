#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           3a_submit_filter_format_true.sh
###############################################################################


set -euo pipefail
source /home1/karatas/ibd_compare/other_scripts/const.sh
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh

# submit formatting for true outputs
log_msg "Detected filter and format specifications
min_bp_len=${MIN_BP_LEN}
min_cm_len=${MIN_CM_LEN}
common_ff_dir=${COMMON_FF_DIR}
chr=${CHR}"
total_tasks=$(( NUM_REPS * ${#TMAPS[@]} ))
sbatch \
    --chdir=/home1/karatas/ibd_compare \
    --job-name=trueFF \
    --array="1-${total_tasks}%${MAX_JOBS}" \
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
