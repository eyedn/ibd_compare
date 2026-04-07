#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           4a_submit_true_vs_hapibd.sh
###############################################################################


set -euo pipefail
source /home1/karatas/ibd_compare/other_scripts/const.sh
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh

total_tasks=$(( ${#TMAPS[@]} * ${#IMAPS[@]} * ${#FILTERING[@]} * ${#POPS[@]} ))

sbatch \
    --chdir=/home1/karatas/ibd_compare \
    --job-name="tvHAP" \
    --array="1-${total_tasks}%${MAX_JOBS}" \
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
total_tasks=${total_tasks}
num_reps=${NUM_REPS}
tvi_dir=${TVI_DIR}

hapibd params:
hap_min_seed=${HAP_MIN_SEED}
hap_max_gap=${HAP_MAX_GAP}
hap_min_extend=${HAP_MIN_EXTEND}
hap_min_output=${HAP_MIN_OUTPUT}
hap_min_markers=${HAP_MIN_MARKERS}
hap_min_mac=${HAP_MIN_MAC}"
