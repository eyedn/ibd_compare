#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           2c_submit_inference_rapid.sh
###############################################################################


set -euo pipefail
source /home1/karatas/ibd_compare/other_scripts/const.sh
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh
source /home1/karatas/ibd_compare/other_scripts/software_profiles.sh
select_software_profile rapid "$@"


# submit ibd inference
log_msg "Detected rapid params
sim_chr=${CHR}
outdir=${OUTDIR}
rapid_map_dir=${RAPID_MAP_DIR}
min_len=${RAP_MIN_LEN} 
window_size=${RAP_WINDOW_SIZE}
num_runs=${RAP_NUM_RUNS}
num_successes=${RAP_NUM_SUCCESSES}"
total_tasks=$(( NUM_REPS * ${#TMAPS[@]} * ${#IMAPS[@]} ))
sbatch \
  --chdir=/home1/karatas/ibd_compare \
  --job-name=rapid \
  --array="1-${total_tasks}%${MAX_JOBS}" \
  --exclude=d11-05,d11-06,d11-07,b23-11,b23-13 \
  job_scripts/inference_rapid.sh \
    "$RAP_MIN_LEN" \
    "$RAP_WINDOW_SIZE" \
    "$RAP_NUM_RUNS" \
    "$RAP_NUM_SUCCESSES" \
    "$VCF_DIR" \
    "$RAPID_MAP_DIR" \
    "$OUTDIR" \
    "$NUM_REPS" \
    "$CHR" \
    "$SIM_ID_CONTRACT_FILE" \
    "$SIM_ID_CONTRACT" \
    "${TMAPS[@]}" \
    -- \
    "${IMAPS[@]}" \
    -- \
    "${POPS[@]}"
log_msg "submitted rapid inference total_tasks=${total_tasks}"
