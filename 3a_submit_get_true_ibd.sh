#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           submit_get_true_ibd.sh
###############################################################################

set -euo pipefail
source /home1/karatas/ibd_compare/other_scripts/const.sh
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh

# run true ibd finding for each tree
log_msg "Detected true ibd specifications
MIN_SPAN=${MIN_SPAN}
MAX_TIME=${MAX_TIME}
OUTDIR=${OUTDIR}
CHR=${CHR}"
TOTAL_TASKS=$(( NUM_REPS * ${#TMAPS[@]} ))
sbatch \
    --chdir=/home1/karatas/ibd_compare \
    --job-name=trueIBD \
    --exclude=d11-05,d11-06,d11-07,b23-11,b23-13 \
    --array="1-${TOTAL_TASKS}%${MAX_JOBS}" \
    "job_scripts/get_true_ibd.sh" \
        "$MIN_SPAN" \
        "$MAX_TIME" \
        "$CHR" \
        "$TREE_DIR" \
        "$TRUE_DIR" \
        "$NUM_REPS" \
        "${TMAPS[@]}" \
        -- \
        "${POPS[@]}"
log_msg "submitted getting true ibd; tasks=${TOTAL_TASKS}"
