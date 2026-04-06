#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           submit_true_vs_inference.sh
###############################################################################

set -euo pipefail
source /home1/karatas/ibd_compare/other_scripts/const.sh
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh

TOTAL_TASKS=$(( ${#TMAPS[@]} * ${#IMAPS[@]} * ${#SOFTWARE[@]} * ${#FILTERING[@]} * ${#POPS[@]} ))

# all software outputs are formatted in this script
# note, specifications of software, filtering, and associated labels are
# specified in the script bellow.
# MAKE SURE CONSISTENT WITH INFERENCE AND FILTERING SCRIPTS
tvi_jobid=$(sbatch \
    --parsable \
    --job-name="tviIBD" \
    --array="1-${TOTAL_TASKS}%${MAX_JOBS}" \
    --exclude=d11-05,d11-06,d11-07,b23-11,b23-13 \
    job_scripts/true_vs_inference.sh \
        "$COMMON_FF_DIR" \
        "$POP_INFO_DIR" \
        "$TVI_DIR" \
        "$CHR" \
        "$SAMPLE_SIZE" \
        "$NUM_REPS" \
        "${TMAPS[@]}" \
        -- \
        "${IMAPS[@]}" \
        -- \
        "${SOFTWARE[@]}" \
        -- \
        "${FILTERING[@]}" \
        -- \
        "${POPS[@]}" )
log_msg "submitted true vs inference (replicates handled inside each task)
tvi_jobid=${tvi_jobid}
TOTAL_TASKS=${TOTAL_TASKS}
NUM_REPS=${NUM_REPS}
TVI_DIR=${TVI_DIR}"

log_msg "submitted combine tvi script; dep on tvi_jobid=${tvi_jobid}"
sbatch \
    --dependency=afterok:${tvi_jobid} \
    --job-name="combTVI" \
    job_scripts/combine_tvi.sh \
        "$TVI_DIR" \
        "${TVI_DIR}/all_stats.tsv" \
        "${TVI_DIR}/all_stats.per_rep.tsv"
