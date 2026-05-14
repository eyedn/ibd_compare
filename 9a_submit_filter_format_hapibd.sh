#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           9a_submit_filter_format_hapibd.sh
###############################################################################


set -euo pipefail
source /home1/karatas/ibd_compare/other_scripts/const.sh
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh

total_tasks=$(( ${#ONEKG_CHRS[@]} * ${#IMAPS[@]} ))

log_msg "Submitting 1kG hapibd filter+format
onekg_ibd_init_inference_dir=${ONEKG_IBD_INIT_INFERENCE_DIR}
onekg_ibd_compare_dir=${ONEKG_IBD_COMPARE_DIR}
onekg_vcf_dir=${ONEKG_VCF_DIR}
populations=${ONEKG_POPS[*]}
chromosomes=${ONEKG_CHRS[*]}
imaps=${IMAPS[*]}
total_tasks=${total_tasks}
cpus_per_task=${NONADMIXED_CPUS_PER_TASK}
mem=${NONADMIXED_MEM}
ibdends_threads=${NONADMIXED_IBDENDS_THREADS}
ibdends_xmx=${NONADMIXED_IBDENDS_XMX}
min_bp_len=${MIN_BP_LEN}
min_cm_len=${MIN_CM_LEN}
hot_bin_size=${HOT_BIN_SIZE}
hot_bin_sd=${HOT_BIN_SD}
min_seed=${HAP_MIN_SEED}
max_gap=${HAP_MAX_GAP}
min_extend=${HAP_MIN_EXTEND}
min_output=${HAP_MIN_OUTPUT}
min_markers=${HAP_MIN_MARKERS}
min_mac=${HAP_MIN_MAC}"

sbatch \
    --chdir=/home1/karatas/ibd_compare \
    --job-name="1kGhapFF" \
    --cpus-per-task="${NONADMIXED_CPUS_PER_TASK}" \
    --mem="${NONADMIXED_MEM}" \
    --export=ALL,IBDENDS_XMX="${NONADMIXED_IBDENDS_XMX}",IBDENDS_THREADS="${NONADMIXED_IBDENDS_THREADS}" \
    --array="1-${total_tasks}%${MAX_JOBS}" \
    --exclude=d11-05,d11-06,d11-07,b23-11,b23-13,b18-08 \
    job_scripts/filter_format_1kG_hapibd.sh \
        "$ONEKG_IBD_INIT_INFERENCE_DIR" \
        "$ONEKG_VCF_DIR" \
        "$ONEKG_IBD_COMPARE_DIR" \
        "$REF_DIR" \
        "$MIN_BP_LEN" \
        "$MIN_CM_LEN" \
        "$HOT_BIN_SIZE" \
        "$HOT_BIN_SD" \
        "$HAP_MIN_SEED" \
        "$HAP_MAX_GAP" \
        "$HAP_MIN_EXTEND" \
        "$HAP_MIN_OUTPUT" \
        "$HAP_MIN_MARKERS" \
        "$HAP_MIN_MAC" \
        "${ONEKG_CHRS[@]}" \
        -- \
        "${IMAPS[@]}" \
        -- \
        "${ONEKG_POPS[@]}"
