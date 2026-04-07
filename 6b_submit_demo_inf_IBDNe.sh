#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           6b_submit_demo_inf_IBDNe.sh
###############################################################################

set -euo pipefail
source /home1/karatas/ibd_compare/other_scripts/const.sh
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh

# run true ibd finding for each tree
log_msg "Detected IBDNe submission specifications
VCF_DIR=${VCF_DIR}
DEMO_DIR=${DEMO_DIR}
CHR=${CHR}
FLAT_RATE=${FLAT_RATE}"
TOTAL_TASKS=$(( NUM_REPS * ${#TMAPS[@]} ))
sbatch \
    --chdir=/home1/karatas/ibd_compare \
    --job-name=demoIBDNe \
    --array="1-${TOTAL_TASKS}%${MAX_JOBS}" \
    "job_scripts/demo_inf_IBDNe.sh" \
        "$VCF_DIR" \
        "$DEMO_DIR" \
        "$CHR" \
        "$FLAT_RATE" \
        "$NUM_REPS" \
        "${TMAPS[@]}" \
        -- \
        "${POPS[@]}"
log_msg "submitted ibdne inference; tasks=${TOTAL_TASKS}"
