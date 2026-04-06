#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           submit_filter_format_hapibd.sh
###############################################################################

set -euo pipefail
source /home1/karatas/ibd_compare/other_scripts/const.sh
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh

log_msg "Submitting hapibd filter+format
MIN_BP_LEN=${MIN_BP_LEN}
MIN_CM_LEN=${MIN_CM_LEN}
HOT_BIN_SIZE=${HOT_BIN_SIZE}
HOT_BIN_SD=${HOT_BIN_SD}

hapibd params: 
HAP_MIN_SEED=${HAP_MIN_SEED} 
HAP_MAX_GAP=${HAP_MAX_GAP} 
HAP_MIN_EXTEND=${HAP_MIN_EXTEND} 
HAP_MIN_OUTPUT=${HAP_MIN_OUTPUT} 
HAP_MIN_MARKERS=${HAP_MIN_MARKERS} 
HAP_MIN_MAC=${HAP_MIN_MAC}
CHR=${CHR}"

TOTAL_TASKS=$(( NUM_REPS * ${#TMAPS[@]} * ${#IMAPS[@]} ))

sbatch \
    --chdir=/home1/karatas/ibd_compare \
    --job-name=hapFF \
    --array="1-${TOTAL_TASKS}%${MAX_JOBS}" \
    /home1/karatas/ibd_compare/job_scripts/filter_format_hapibd.sh \
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
        "$HAP_MIN_SEED" \
        "$HAP_MAX_GAP" \
        "$HAP_MIN_EXTEND" \
        "$HAP_MIN_OUTPUT" \
        "$HAP_MIN_MARKERS" \
        "$HAP_MIN_MAC" \
        "${TMAPS[@]}" \
        -- \
        "${IMAPS[@]}" \
        -- \
        "${POPS[@]}"

log_msg "submitted hapibd filter+format jobs; TOTAL_TASKS=${TOTAL_TASKS}"
