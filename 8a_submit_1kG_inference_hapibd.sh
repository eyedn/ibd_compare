#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           8a_submit_1kG_inference_hapibd.sh
###############################################################################


set -euo pipefail
source /home1/karatas/ibd_compare/other_scripts/const.sh
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh

total_tasks=$(( ${#ONEKG_CHRS[@]} * ${#IMAPS[@]} ))

log_msg "Submitting 1kG hapibd inference
onekg_ibd_compare_dir=${ONEKG_IBD_COMPARE_DIR}
onekg_vcf_dir=${ONEKG_VCF_DIR}
populations=${ONEKG_POPS[*]}
chromosomes=${ONEKG_CHRS[*]}
imaps=${IMAPS[*]}
total_tasks=${total_tasks}
min_seed=${HAP_MIN_SEED}
max_gap=${HAP_MAX_GAP}
min_extend=${HAP_MIN_EXTEND}
min_output=${HAP_MIN_OUTPUT}
min_markers=${HAP_MIN_MARKERS}
min_mac=${HAP_MIN_MAC}"

sbatch \
    --chdir=/home1/karatas/ibd_compare \
    --job-name="1kGhapibd" \
    --array="1-${total_tasks}%${MAX_JOBS}" \
    --exclude=d11-05,d11-06,d11-07,b23-11,b23-13 \
    job_scripts/inference_1kG_hapibd.sh \
        "$HAP_MIN_SEED" \
        "$HAP_MAX_GAP" \
        "$HAP_MIN_EXTEND" \
        "$HAP_MIN_OUTPUT" \
        "$HAP_MIN_MARKERS" \
        "$HAP_MIN_MAC" \
        "$ONEKG_IBD_COMPARE_DIR" \
        "$ONEKG_VCF_DIR" \
        "$REF_DIR" \
        "${ONEKG_CHRS[@]}" \
        -- \
        "${IMAPS[@]}" \
        -- \
        "${ONEKG_POPS[@]}"
