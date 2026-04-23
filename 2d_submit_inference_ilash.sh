#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           2d_submit_inference_ilash.sh
###############################################################################


set -euo pipefail
source /home1/karatas/ibd_compare/other_scripts/const.sh
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh

log_msg "Detected ilash params
sim_chr=${CHR}
outdir=${OUTDIR}
ped_dir=${PED_DIR}
ilash_map_dir=${ILASH_MAP_DIR}
auto_slice=${ILASH_AUTO_SLICE}
slice_size=${ILASH_SLICE_SIZE}
step_size=${ILASH_STEP_SIZE}
perm_count=${ILASH_PERM_COUNT}
shingle_size=${ILASH_SHINGLE_SIZE}
shingle_overlap=${ILASH_SHINGLE_OVERLAP}
bucket_count=${ILASH_BUCKET_COUNT}
match_threshold=${ILASH_MATCH_THRESHOLD}
interest_threshold=${ILASH_INTEREST_THRESHOLD}
min_length=${ILASH_MIN_LENGTH}
slice_length=${ILASH_SLICE_LENGTH}
cm_overlap=${ILASH_CM_OVERLAP}
minhash_threshold=${ILASH_MINHASH_THRESHOLD}"

total_tasks=$(( NUM_REPS * ${#TMAPS[@]} * ${#IMAPS[@]} ))

sbatch \
    --chdir=/home1/karatas/ibd_compare \
    --job-name=ilash \
    --array="1-${total_tasks}%${MAX_JOBS}" \
    --exclude=d11-05,d11-06,d11-07,b23-11,b23-13 \
    job_scripts/inference_ilash.sh \
        "$ILASH_AUTO_SLICE" \
        "$ILASH_SLICE_SIZE" \
        "$ILASH_STEP_SIZE" \
        "$ILASH_PERM_COUNT" \
        "$ILASH_SHINGLE_SIZE" \
        "$ILASH_SHINGLE_OVERLAP" \
        "$ILASH_BUCKET_COUNT" \
        "$ILASH_MATCH_THRESHOLD" \
        "$ILASH_INTEREST_THRESHOLD" \
        "$ILASH_MIN_LENGTH" \
        "$ILASH_SLICE_LENGTH" \
        "$ILASH_CM_OVERLAP" \
        "$ILASH_MINHASH_THRESHOLD" \
        "$OUTDIR" \
        "$PED_DIR" \
        "$ILASH_MAP_DIR" \
        "$NUM_REPS" \
        "$CHR" \
        "${TMAPS[@]}" \
        -- \
        "${IMAPS[@]}" \
        -- \
        "${POPS[@]}"

log_msg "submitted ilash inference total_tasks=${total_tasks}"
