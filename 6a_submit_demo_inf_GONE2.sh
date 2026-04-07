#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           6a_submit_demo_inf_GONE2.sh
###############################################################################


set -euo pipefail
source /home1/karatas/ibd_compare/other_scripts/const.sh
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh

# run true ibd finding for each tree
log_msg "Detected GONE2 submission specifications
vcf_dir=${VCF_DIR}
demo_dir=${DEMO_DIR}
chr=${CHR}
flat_rate=${FLAT_RATE}"
total_tasks=$(( NUM_REPS * ${#TMAPS[@]} ))
sbatch \
    --chdir=/home1/karatas/ibd_compare \
    --job-name=demoGONE2 \
    --array="1-${total_tasks}%${MAX_JOBS}" \
    "job_scripts/demo_inf_GONE2.sh" \
        "$VCF_DIR" \
        "$DEMO_DIR" \
        "$CHR" \
        "$FLAT_RATE" \
        "$NUM_REPS" \
        "${TMAPS[@]}" \
        -- \
        "${POPS[@]}"
log_msg "submitted gone2 inference; tasks=${total_tasks}"
