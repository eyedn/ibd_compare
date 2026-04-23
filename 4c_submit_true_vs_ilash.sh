#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           4c_submit_true_vs_ilash.sh
###############################################################################


set -euo pipefail
source /home1/karatas/ibd_compare/other_scripts/const.sh
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh

total_tasks=$(( ${#TMAPS[@]} * ${#IMAPS[@]} * ${#FILTERING[@]} * ${#POPS[@]} ))

sbatch \
    --chdir=/home1/karatas/ibd_compare \
    --job-name="tvILA" \
    --array="1-${total_tasks}%${MAX_JOBS}" \
    --exclude=d11-05,d11-06,d11-07,b23-11,b23-13 \
    /home1/karatas/ibd_compare/job_scripts/true_vs_ilash.sh \
        "$COMMON_FF_DIR" \
        "$POP_INFO_DIR" \
        "$TVI_DIR" \
        "$CHR" \
        "$SAMPLE_SIZE" \
        "$NUM_REPS" \
        "$ILASH_AUTO_SLICE" \
        "$ILASH_MIN_LENGTH" \
        "$ILASH_SLICE_LENGTH" \
        "$ILASH_CM_OVERLAP" \
        "${TMAPS[@]}" \
        -- \
        "${IMAPS[@]}" \
        -- \
        "${FILTERING[@]}" \
        -- \
        "${POPS[@]}"

log_msg "submitted ilash true-vs-inference
total_tasks=${total_tasks}
num_reps=${NUM_REPS}
tvi_dir=${TVI_DIR}

ilash params:
auto_slice=${ILASH_AUTO_SLICE}
min_length=${ILASH_MIN_LENGTH}
slice_length=${ILASH_SLICE_LENGTH}
cm_overlap=${ILASH_CM_OVERLAP}"
