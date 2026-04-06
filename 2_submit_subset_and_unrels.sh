#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           submit_subset_and_unrels.sh
###############################################################################

set -euo pipefail
source /home1/karatas/ibd_compare/other_scripts/const.sh
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh

# subset vcf file by population
TOTAL_TASKS=$(( NUM_REPS * ${#TMAPS[@]} ))
subset_jobid=$(sbatch \
    --parsable \
    --chdir=/home1/karatas/ibd_compare \
    --exclude=d11-05,d11-06,d11-07,b23-11,b23-13,b18-08 \
    --job-name="sub" \
    --array="1-${TOTAL_TASKS}%${MAX_JOBS}" \
    "job_scripts/subset_ts.sh" \
        "$TREE_DIR" \
        "$VCF_DIR" \
        "$CHR" \
        "$NUM_REPS" \
        "${TMAPS[@]}" \
        -- \
        "${POPS[@]}"
)
log_msg "subset ts array for each pop and map and create vcfs; jid=${subset_jobid}
TOTAL_TASKS=${TOTAL_TASKS}"

# get unrelated individuals from subsetted vcf files
unrel_jobid=$(sbatch \
    --parsable \
    --chdir=/home1/karatas/ibd_compare \
    --exclude=d11-05,d11-06,d11-07,b23-11,b23-13,b18-08 \
    --dependency=afterok:${subset_jobid} \
    --job-name="unrel" \
    --array="1-${TOTAL_TASKS}%${MAX_JOBS}" \
    "job_scripts/get_unrels.sh" \
        "$KIN_CUTOFF" \
        "$VCF_DIR" \
        "$POP_INFO_DIR" \
        "$NUM_REPS" \
        "${TMAPS[@]}" \
        -- \
        "${POPS[@]}"
)
log_msg "subset vcf unrels for each pop and map; jid=${unrel_jobid}; dep. on ${subset_jobid}
TOTAL_TASKS=${TOTAL_TASKS}"
