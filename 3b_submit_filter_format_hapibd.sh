#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           3b_submit_filter_format_hapibd.sh
###############################################################################


set -euo pipefail
source /home1/karatas/ibd_compare/other_scripts/const.sh
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh

log_msg "Submitting hapibd filter+format
array dimensions: tmap x imap x seed x pop (one pop per task)
min_bp_len=${MIN_BP_LEN}
min_cm_len=${MIN_CM_LEN}
hot_bin_size=${HOT_BIN_SIZE}
hot_bin_sd=${HOT_BIN_SD}

hapibd params: 
hap_min_seed=${HAP_MIN_SEED} 
hap_max_gap=${HAP_MAX_GAP} 
hap_min_extend=${HAP_MIN_EXTEND} 
hap_min_output=${HAP_MIN_OUTPUT} 
hap_min_markers=${HAP_MIN_MARKERS} 
hap_min_mac=${HAP_MIN_MAC}
chr=${CHR}"

total_tasks=$(( NUM_REPS * ${#TMAPS[@]} * ${#IMAPS[@]} * ${#POPS[@]} ))

sbatch \
    --chdir=/home1/karatas/ibd_compare \
    --job-name=hapFF \
    --array="1-${total_tasks}%${MAX_JOBS}" \
    /home1/karatas/ibd_compare/job_scripts/filter_format_hapibd.sh \
        "$OUTDIR" \
        "$VCF_DIR" \
        "$COMMON_FF_DIR" \
        "$REF_DIR" \
        "$CHR" \
        "$MIN_BP_LEN" \
        "$MIN_CM_LEN" \
        "$HOT_BIN_SIZE" \
        "$HOT_BIN_SD" \
        "$NUM_REPS" \
        "$HAP_MIN_SEED" \
        "$HAP_MAX_GAP" \
        "$HAP_MIN_EXTEND" \
        "$HAP_MIN_OUTPUT" \
        "$HAP_MIN_MARKERS" \
        "$HAP_MIN_MAC" \
        "${TMAPS[@]}" \
        -- \
        "${IMAPS[@]}" \
        -- \
        "${POPS[@]}"

log_msg "submitted hapibd filter+format jobs; total_tasks=${total_tasks}"
