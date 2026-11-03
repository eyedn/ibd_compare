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
source /home1/karatas/ibd_compare/other_scripts/software_profiles.sh
select_software_profile hapibd "$@"


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
ibdends_min_maf=${IBDENDS_MIN_MAF}
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

# eachs path has its own MAX_JOBS cap; the two roots run independently.
hot_jid=""
ibdends_jid=""
for filter_path in hot ibdends hotibdends ibdendshot; do
    dependency_args=()
    dependency="none"
    case "$filter_path" in
        hotibdends) dependency="afterok:${hot_jid}" ;;
        ibdendshot) dependency="afterok:${ibdends_jid}" ;;
    esac
    if [[ "$dependency" != none ]]; then
        dependency_args=("--dependency=${dependency}")
    fi
    job_id=$(sbatch \
        --parsable \
        ${dependency_args[@]+"${dependency_args[@]}"} \
        --chdir=/home1/karatas/ibd_compare \
        --job-name="1kGhapFF_${filter_path}" \
        --cpus-per-task="${NONADMIXED_CPUS_PER_TASK}" \
        --mem="${NONADMIXED_MEM}" \
        --export=ALL,IBDENDS_XMX="${NONADMIXED_IBDENDS_XMX}",IBDENDS_THREADS="${NONADMIXED_IBDENDS_THREADS}",IBDENDS_MIN_MAF="${IBDENDS_MIN_MAF}" \
        --array="1-${total_tasks}%${MAX_JOBS}" \
        --exclude=d11-05,d11-06,d11-07,b23-11,b23-13,b18-08 \
        job_scripts/filter_format_1kG_hapibd.sh \
            "$filter_path" \
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
    )
    job_id="${job_id%%;*}"
    [[ "$job_id" =~ ^[0-9]+$ ]] || {
        echo "ERROR: invalid sbatch job ID: $job_id" >&2
        exit 1
    }
    case "$filter_path" in
        hot) hot_jid="$job_id" ;;
        ibdends) ibdends_jid="$job_id" ;;
        hotibdends) hotibdends_jid="$job_id" ;;
        ibdendshot) ibdendshot_jid="$job_id" ;;
    esac
    log_msg "submitted 1kG hapibd path=$filter_path; jid=$job_id; dep=$dependency"
done
