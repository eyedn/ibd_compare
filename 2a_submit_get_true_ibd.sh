#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           2a_submit_get_true_ibd.sh
###############################################################################


set -euo pipefail
source /home1/karatas/ibd_compare/other_scripts/const.sh
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh

# run true ibd finding for each tree
log_msg "Detected true ibd specifications
min_span=${MIN_SPAN}
max_time=${MAX_TIME}
outdir=${OUTDIR}
chr=${CHR}"
total_tasks=$(( NUM_REPS * ${#TMAPS[@]} ))
sbatch \
    --chdir=/home1/karatas/ibd_compare \
    --job-name=trueIBD \
    --exclude=d11-05,d11-06,d11-07,b23-11,b23-13 \
    --array="1-${total_tasks}%${MAX_JOBS}" \
    "job_scripts/get_true_ibd.sh" \
        "$MIN_SPAN" \
        "$MAX_TIME" \
        "$CHR" \
        "$OOA_TREE_DIR" \
        "$TRUE_DIR" \
        "$SAMPLE_SIZE" \
        "$ANALYSIS_SAMPLE_SIZE" \
        "$NUM_REPS" \
        "$SIM_ID_CONTRACT_FILE" \
        "$SIM_ID_CONTRACT" \
        "${TMAPS[@]}" \
        -- \
        "${POPS[@]}"
log_msg "submitted getting true ibd; tasks=${total_tasks}"
