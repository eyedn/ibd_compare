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

standard_pops=()
admix_pops=()
for pop in "${POPS[@]}"; do
    if [[ "$pop" == "ADMIX" ]]; then
        admix_pops+=( "$pop" )
    else
        standard_pops+=( "$pop" )
    fi
done

if [[ ${#standard_pops[@]} -gt 0 ]]; then
    total_tasks=$(( NUM_REPS * ${#TMAPS[@]} * ${#IMAPS[@]} * ${#standard_pops[@]} ))
    log_msg "Submitting rapid standard-tier jobs
populations=${standard_pops[*]}
cpus_per_task=4
mem=32G
ibdends_threads=4
ibdends_xmx=20G
total_tasks=${total_tasks}"

    sbatch \
        --chdir=/home1/karatas/ibd_compare \
        --job-name=rapFFstd \
        --cpus-per-task=4 \
        --mem=32G \
        --export=ALL,IBDENDS_XMX=20G,IBDENDS_THREADS=4 \
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
            "${standard_pops[@]}"
fi

if [[ ${#admix_pops[@]} -gt 0 ]]; then
    total_tasks=$(( NUM_REPS * ${#TMAPS[@]} * ${#IMAPS[@]} * ${#admix_pops[@]} ))
    log_msg "Submitting rapid ADMIX-tier jobs
populations=${admix_pops[*]}
cpus_per_task=4
mem=128G
ibdends_threads=4
ibdends_xmx=96G
total_tasks=${total_tasks}"

    sbatch \
        --chdir=/home1/karatas/ibd_compare \
        --job-name=rapFFadx \
        --cpus-per-task=4 \
        --mem=128G \
        --export=ALL,IBDENDS_XMX=96G,IBDENDS_THREADS=4 \
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
            "${admix_pops[@]}"
fi
