#!/usr/bin/env bash

#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           submit_true_vs_hapibd.sh
###############################################################################

set -euo pipefail
source /home1/karatas/ibd_compare/other_scripts/const.sh
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh

TOTAL_TASKS=$(( ${#TMAPS[@]} * ${#IMAPS[@]} * ${#FILTERING[@]} * ${#POPS[@]} ))

sbatch \
    --chdir=/home1/karatas/ibd_compare \
    --job-name="tvHAP" \
    --array="1-${TOTAL_TASKS}%${MAX_JOBS}" \
    --exclude=d11-05,d11-06,d11-07,b23-11,b23-13 \
    /home1/karatas/ibd_compare/job_scripts/true_vs_hapibd.sh \
        "$COMMON_FF_DIR" \
        "$POP_INFO_DIR" \
        "$TVI_DIR" \
        "$CHR" \
        "$SAMPLE_SIZE" \
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
        "${FILTERING[@]}" \
        -- \
    "${POPS[@]}"
log_msg "submitted hapibd true-vs-inference
TOTAL_TASKS=${TOTAL_TASKS}
NUM_REPS=${NUM_REPS}
TVI_DIR=${TVI_DIR}

hapibd params:
HAP_MIN_SEED=${HAP_MIN_SEED}
HAP_MAX_GAP=${HAP_MAX_GAP}
HAP_MIN_EXTEND=${HAP_MIN_EXTEND}
HAP_MIN_OUTPUT=${HAP_MIN_OUTPUT}
HAP_MIN_MARKERS=${HAP_MIN_MARKERS}
HAP_MIN_MAC=${HAP_MIN_MAC}"
