#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           submit_sim.sh
###############################################################################

set -euo pipefail
source /home1/karatas/ibd_compare/other_scripts/const.sh
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh

# NOTE:  sim. specs in /home1/karatas/ibd_compare/const.sh
TOTAL_TASKS=$(( NUM_REPS * ${#TMAPS[@]} ))
log_msg "Detected sim specifications
OUTDIR=${OUTDIR}
NUM_REPS=${NUM_REPS}
SAMPLE_SIZE=${SAMPLE_SIZE}
TREE_DIR=${TREE_DIR}
MSPRIME_MODEL=${MSPRIME_MODEL}
DEMO_MODEL=${DEMO_MODEL}
CHR=${CHR}
TOTAL_TASKS=${TOTAL_TASKS} (per job array)"

# submit sims
TOTAL_TASKS=$(( NUM_REPS * ${#TMAPS[@]} ))
sim_jid=$(sbatch \
    --parsable \
    --chdir=/home1/karatas/ibd_compare \
    --exclude=d11-05,d11-06,d11-07,b23-11,b23-13,b18-08 \
    --dependency=afterok:${subset_jid} \
    --job-name="sim" \
    --array="1-${TOTAL_TASKS}%${MAX_JOBS}" \
    "job_scripts/sim_ooaa.sh" \
        "$TREE_DIR" "$SAMPLE_SIZE" "$NUM_REPS" \
        "$MSPRIME_MODEL" "$DEMO_MODEL" "$CHR" \
        "${TMAPS[@]}"
)
log_msg "submitted combined sim array; jid=${sim_jid}"

# subset vcf file by population
TOTAL_TASKS=$(( NUM_REPS * ${#TMAPS[@]} ))
subset_jid=$(sbatch \
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
log_msg "subset ts array for each pop and create vcfs; jid=${subset_jid}; dep on ${sim_jid}"

# get unrelated individuals from subsetted vcf files
unrel_jid=$(sbatch \
    --parsable \
    --chdir=/home1/karatas/ibd_compare \
    --exclude=d11-05,d11-06,d11-07,b23-11,b23-13,b18-08 \
    --dependency=afterok:${subset_jid} \
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
log_msg "subset vcf unrels for each pop and map; jid=${unrel_jid}; dep. on ${subset_jid}"
