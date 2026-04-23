#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           3c_submit_filter_format_rapid.sh
###############################################################################


set -euo pipefail
source /home1/karatas/ibd_compare/other_scripts/const.sh
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh

log_msg "Submitting rapid filter+format
array dimensions: tmap x imap x seed x pop (one pop per task)
min_bp_len=${MIN_BP_LEN}
min_cm_len=${MIN_CM_LEN}
hot_bin_size=${HOT_BIN_SIZE}
hot_bin_sd=${HOT_BIN_SD}

rapid params: 
rap_min_len=${RAP_MIN_LEN} 
rap_window_size=${RAP_WINDOW_SIZE} 
rap_num_runs=${RAP_NUM_RUNS} 
rap_num_successes=${RAP_NUM_SUCCESSES}
rapid_map_dir=${RAPID_MAP_DIR}
chr=${CHR}"

total_tasks=$(( NUM_REPS * ${#TMAPS[@]} * ${#IMAPS[@]} * ${#NONADMIXED_POP_LABELS[@]} ))
log_msg "Submitting rapid nonadmixed-tier jobs
populations=${NONADMIXED_POP_LABELS[*]}
cpus_per_task=${NONADMIXED_CPUS_PER_TASK}
mem=${NONADMIXED_MEM}
ibdends_threads=${NONADMIXED_IBDENDS_THREADS}
ibdends_xmx=${NONADMIXED_IBDENDS_XMX}
total_tasks=${total_tasks}"

sbatch \
    --chdir=/home1/karatas/ibd_compare \
    --job-name=rapFFnAdx \
    --cpus-per-task="${NONADMIXED_CPUS_PER_TASK}" \
    --mem="${NONADMIXED_MEM}" \
    --export=ALL,IBDENDS_XMX="${NONADMIXED_IBDENDS_XMX}",IBDENDS_THREADS="${NONADMIXED_IBDENDS_THREADS}" \
    --array="1-${total_tasks}%${MAX_JOBS}" \
    /home1/karatas/ibd_compare/job_scripts/filter_format_rapid.sh \
        "$OUTDIR" \
        "$VCF_DIR" \
        "$RAPID_MAP_DIR" \
        "$COMMON_FF_DIR" \
        "$REF_DIR" \
        "$CHR" \
        "$MIN_BP_LEN" \
        "$MIN_CM_LEN" \
        "$HOT_BIN_SIZE" \
        "$HOT_BIN_SD" \
        "$NUM_REPS" \
        "$RAP_MIN_LEN" \
        "$RAP_WINDOW_SIZE" \
        "$RAP_NUM_RUNS" \
        "$RAP_NUM_SUCCESSES" \
        "${TMAPS[@]}" \
        -- \
        "${IMAPS[@]}" \
        -- \
        "${NONADMIXED_POP_LABELS[@]}"

total_tasks=$(( NUM_REPS * ${#TMAPS[@]} * ${#IMAPS[@]} ))
log_msg "Submitting rapid ADMIX-tier jobs
populations=${ADMIXED_POP_LABEL}
cpus_per_task=${ADMIXED_CPUS_PER_TASK}
mem=${ADMIXED_MEM}
ibdends_threads=${ADMIXED_IBDENDS_THREADS}
ibdends_xmx=${ADMIXED_IBDENDS_XMX}
total_tasks=${total_tasks}"

sbatch \
    --chdir=/home1/karatas/ibd_compare \
    --job-name=rapFFAdx \
    --cpus-per-task="${ADMIXED_CPUS_PER_TASK}" \
    --mem="${ADMIXED_MEM}" \
    --export=ALL,IBDENDS_XMX="${ADMIXED_IBDENDS_XMX}",IBDENDS_THREADS="${ADMIXED_IBDENDS_THREADS}" \
    --array="1-${total_tasks}%${MAX_JOBS}" \
    /home1/karatas/ibd_compare/job_scripts/filter_format_rapid.sh \
        "$OUTDIR" \
        "$VCF_DIR" \
        "$RAPID_MAP_DIR" \
        "$COMMON_FF_DIR" \
        "$REF_DIR" \
        "$CHR" \
        "$MIN_BP_LEN" \
        "$MIN_CM_LEN" \
        "$HOT_BIN_SIZE" \
        "$HOT_BIN_SD" \
        "$NUM_REPS" \
        "$RAP_MIN_LEN" \
        "$RAP_WINDOW_SIZE" \
        "$RAP_NUM_RUNS" \
        "$RAP_NUM_SUCCESSES" \
        "${TMAPS[@]}" \
        -- \
        "${IMAPS[@]}" \
        -- \
        "$ADMIXED_POP_LABEL"
