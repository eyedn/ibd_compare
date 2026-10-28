#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           4b_submit_true_vs_rapid.sh
###############################################################################


set -euo pipefail
source /home1/karatas/ibd_compare/other_scripts/const.sh
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh
source /home1/karatas/ibd_compare/other_scripts/software_profiles.sh
select_software_profile rapid "$@"


total_tasks=$(( ${#TMAPS[@]} * ${#IMAPS[@]} * ${#FILTERING[@]} * ${#POPS[@]} ))

sbatch \
    --chdir=/home1/karatas/ibd_compare \
    --job-name="tvRAP" \
    --array="1-${total_tasks}%${MAX_JOBS}" \
    --exclude=d11-05,d11-06,d11-07,b23-11,b23-13 \
    /home1/karatas/ibd_compare/job_scripts/true_vs_rapid.sh \
        "$COMMON_FF_DIR" \
        "$POP_INFO_DIR" \
        "$TVI_DIR" \
        "$CHR" \
        "$ANALYSIS_SAMPLE_SIZE" \
        "$NUM_REPS" \
        "$RAP_MIN_LEN" \
        "$RAP_WINDOW_SIZE" \
        "$RAP_NUM_RUNS" \
        "$RAP_NUM_SUCCESSES" \
        "$SIM_ID_CONTRACT_FILE" \
        "$SIM_ID_CONTRACT" \
        "${TMAPS[@]}" \
        -- \
        "${IMAPS[@]}" \
        -- \
        "${FILTERING[@]}" \
        -- \
        "${POPS[@]}"
log_msg "submitted rapid true-vs-inference
total_tasks=${total_tasks}
num_reps=${NUM_REPS}
tvi_dir=${TVI_DIR}

rapid params:
rap_min_len=${RAP_MIN_LEN}
rap_window_size=${RAP_WINDOW_SIZE}
rap_num_runs=${RAP_NUM_RUNS}
rap_num_successe=${RAP_NUM_SUCCESSES}"
