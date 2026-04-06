#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           submit_filter_format_rapid.sh
###############################################################################

set -euo pipefail
source /home1/karatas/ibd_compare/other_scripts/const.sh
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh

log_msg "Submitting rapid filter+format
MIN_BP_LEN=${MIN_BP_LEN}
MIN_CM_LEN=${MIN_CM_LEN}
HOT_BIN_SIZE=${HOT_BIN_SIZE}
HOT_BIN_SD=${HOT_BIN_SD}

rapid params: 
RAP_MIN_LEN=${RAP_MIN_LEN} 
RAP_WINDOW_SIZE=${RAP_WINDOW_SIZE} 
RAP_NUM_RUNS=${RAP_NUM_RUNS} 
RAP_NUM_SUCCESSES=${RAP_NUM_SUCCESSES}
CHR=${CHR}"

TOTAL_TASKS=$(( NUM_REPS * ${#TMAPS[@]} * ${#IMAPS[@]} ))

sbatch \
    --chdir=/home1/karatas/ibd_compare \
    --job-name=rapFF \
    --array="1-${TOTAL_TASKS}%${MAX_JOBS}" \
    /home1/karatas/ibd_compare/job_scripts/filter_format_rapid.sh \
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
        "$RAP_MIN_LEN" \
        "$RAP_WINDOW_SIZE" \
        "$RAP_NUM_RUNS" \
        "$RAP_NUM_SUCCESSES" \
        "${TMAPS[@]}" \
        -- \
        "${IMAPS[@]}" \
        -- \
        "${POPS[@]}"

log_msg "submitted rapid filter+format jobs; TOTAL_TASKS=${TOTAL_TASKS}"
