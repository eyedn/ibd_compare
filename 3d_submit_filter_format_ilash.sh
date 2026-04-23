#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           3d_submit_filter_format_ilash.sh
###############################################################################


set -euo pipefail
source /home1/karatas/ibd_compare/other_scripts/const.sh
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh

log_msg "Submitting ilash filter+format
array dimensions: tmap x imap x seed x pop (one pop per task)
min_bp_len=${MIN_BP_LEN}
min_cm_len=${MIN_CM_LEN}
hot_bin_size=${HOT_BIN_SIZE}
hot_bin_sd=${HOT_BIN_SD}

ilash params:
auto_slice=${ILASH_AUTO_SLICE}
min_length=${ILASH_MIN_LENGTH}
slice_length=${ILASH_SLICE_LENGTH}
cm_overlap=${ILASH_CM_OVERLAP}
chr=${CHR}"

total_tasks=$(( NUM_REPS * ${#TMAPS[@]} * ${#IMAPS[@]} * ${#NONADMIXED_POP_LABELS[@]} ))
log_msg "Submitting ilash nonadmixed-tier jobs
populations=${NONADMIXED_POP_LABELS[*]}
cpus_per_task=${NONADMIXED_CPUS_PER_TASK}
mem=${NONADMIXED_MEM}
ibdends_threads=${NONADMIXED_IBDENDS_THREADS}
ibdends_xmx=${NONADMIXED_IBDENDS_XMX}
total_tasks=${total_tasks}"

sbatch \
    --chdir=/home1/karatas/ibd_compare \
    --job-name=ilaFFnAdx \
    --cpus-per-task="${NONADMIXED_CPUS_PER_TASK}" \
    --mem="${NONADMIXED_MEM}" \
    --export=ALL,IBDENDS_XMX="${NONADMIXED_IBDENDS_XMX}",IBDENDS_THREADS="${NONADMIXED_IBDENDS_THREADS}" \
    --array="1-${total_tasks}%${MAX_JOBS}" \
    /home1/karatas/ibd_compare/job_scripts/filter_format_ilash.sh \
        "$OUTDIR" \
        "$VCF_DIR" \
        "$ILASH_MAP_DIR" \
        "$COMMON_FF_DIR" \
        "$CHR" \
        "$MIN_BP_LEN" \
        "$MIN_CM_LEN" \
        "$HOT_BIN_SIZE" \
        "$HOT_BIN_SD" \
        "$NUM_REPS" \
        "$ILASH_AUTO_SLICE" \
        "$ILASH_MIN_LENGTH" \
        "$ILASH_SLICE_LENGTH" \
        "$ILASH_CM_OVERLAP" \
        "${TMAPS[@]}" \
        -- \
        "${IMAPS[@]}" \
        -- \
        "${NONADMIXED_POP_LABELS[@]}"

total_tasks=$(( NUM_REPS * ${#TMAPS[@]} * ${#IMAPS[@]} ))
log_msg "Submitting ilash ADMIX-tier jobs
populations=${ADMIXED_POP_LABEL}
cpus_per_task=${ADMIXED_CPUS_PER_TASK}
mem=${ADMIXED_MEM}
ibdends_threads=${ADMIXED_IBDENDS_THREADS}
ibdends_xmx=${ADMIXED_IBDENDS_XMX}
total_tasks=${total_tasks}"

sbatch \
    --chdir=/home1/karatas/ibd_compare \
    --job-name=ilaFFAdx \
    --cpus-per-task="${ADMIXED_CPUS_PER_TASK}" \
    --mem="${ADMIXED_MEM}" \
    --export=ALL,IBDENDS_XMX="${ADMIXED_IBDENDS_XMX}",IBDENDS_THREADS="${ADMIXED_IBDENDS_THREADS}" \
    --array="1-${total_tasks}%${MAX_JOBS}" \
    /home1/karatas/ibd_compare/job_scripts/filter_format_ilash.sh \
        "$OUTDIR" \
        "$VCF_DIR" \
        "$ILASH_MAP_DIR" \
        "$COMMON_FF_DIR" \
        "$CHR" \
        "$MIN_BP_LEN" \
        "$MIN_CM_LEN" \
        "$HOT_BIN_SIZE" \
        "$HOT_BIN_SD" \
        "$NUM_REPS" \
        "$ILASH_AUTO_SLICE" \
        "$ILASH_MIN_LENGTH" \
        "$ILASH_SLICE_LENGTH" \
        "$ILASH_CM_OVERLAP" \
        "${TMAPS[@]}" \
        -- \
        "${IMAPS[@]}" \
        -- \
        "$ADMIXED_POP_LABEL"
