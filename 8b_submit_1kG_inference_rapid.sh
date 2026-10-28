#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           8b_submit_1kG_inference_rapid.sh
###############################################################################


set -euo pipefail
source /home1/karatas/ibd_compare/other_scripts/const.sh
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh
source /home1/karatas/ibd_compare/other_scripts/software_profiles.sh
select_software_profile rapid "$@"


total_tasks=$(( ${#ONEKG_CHRS[@]} * ${#IMAPS[@]} ))

log_msg "Submitting 1kG rapid inference
onekg_ibd_init_inference_dir=${ONEKG_IBD_INIT_INFERENCE_DIR}
onekg_vcf_dir=${ONEKG_VCF_DIR}
onekg_rapid_map_dir=${ONEKG_RAPID_MAP_DIR}
populations=${ONEKG_POPS[*]}
chromosomes=${ONEKG_CHRS[*]}
imaps=${IMAPS[*]}
total_tasks=${total_tasks}
min_len=${RAP_MIN_LEN}
window_size=${RAP_WINDOW_SIZE}
num_runs=${RAP_NUM_RUNS}
num_successes=${RAP_NUM_SUCCESSES}"

sbatch \
    --chdir=/home1/karatas/ibd_compare \
    --job-name="1kGrapid" \
    --array="1-${total_tasks}%${MAX_JOBS}" \
    --exclude=d11-05,d11-06,d11-07,b23-11,b23-13 \
    job_scripts/inference_1kG_rapid.sh \
        "$RAP_MIN_LEN" \
        "$RAP_WINDOW_SIZE" \
        "$RAP_NUM_RUNS" \
        "$RAP_NUM_SUCCESSES" \
        "$ONEKG_VCF_DIR" \
        "$ONEKG_RAPID_MAP_DIR" \
        "$ONEKG_IBD_INIT_INFERENCE_DIR" \
        "${ONEKG_CHRS[@]}" \
        -- \
        "${IMAPS[@]}" \
        -- \
        "${ONEKG_POPS[@]}"
