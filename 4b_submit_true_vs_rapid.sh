#!/usr/bin/env bash

#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           submit_true_vs_rapid.sh
###############################################################################

set -euo pipefail
source /home1/karatas/ibd_compare/other_scripts/const.sh
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh

TOTAL_TASKS=$(( ${#TMAPS[@]} * ${#IMAPS[@]} * ${#FILTERING[@]} * ${#POPS[@]} ))

sbatch \
    --chdir=/home1/karatas/ibd_compare \
    --job-name="tvRAP" \
    --array="1-${TOTAL_TASKS}%${MAX_JOBS}" \
    --exclude=d11-05,d11-06,d11-07,b23-11,b23-13 \
    /home1/karatas/ibd_compare/job_scripts/true_vs_rapid.sh \
        "$COMMON_FF_DIR" \
        "$POP_INFO_DIR" \
        "$TVI_DIR" \
        "$CHR" \
        "$SAMPLE_SIZE" \
        "$NUM_REPS" \
        "$RAP_MIN_LEN" \
        "$RAP_WINDOW_SIZE" \
        "$RAP_NUM_RUNS" \
        "$RAP_NUM_SUCCESSES" \
        "${TMAPS[@]}" \
        -- \
        "${IMAPS[@]}" \
        -- \
        "${FILTERING[@]}" \
        -- \
        "${POPS[@]}"
log_msg "submitted rapid true-vs-inference
TOTAL_TASKS=${TOTAL_TASKS}
NUM_REPS=${NUM_REPS}
TVI_DIR=${TVI_DIR}

rapid params:
RAP_MIN_LEN=${RAP_MIN_LEN}
RAP_WINDOW_SIZE=${RAP_WINDOW_SIZE}
RAP_NUM_RUNS=${RAP_NUM_RUNS}
RAP_NUM_SUCCESSE=${RAP_NUM_SUCCESSES}"
