#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           7_submit_1kG_prep.sh
###############################################################################


set -euo pipefail
source /home1/karatas/ibd_compare/other_scripts/const.sh
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh

# shared array sizes
vcf_tasks="${#ONEKG_CHRS[@]}"
rapid_map_tasks=$(( ${#ONEKG_CHRS[@]} * ${#IMAPS[@]} * ${#ONEKG_POPS[@]} ))

log_msg "Submitting 1kG prep
onekg_ibd_compare_dir=${ONEKG_IBD_COMPARE_DIR}
onekg_source_vcf_dir=${ONEKG_SOURCE_VCF_DIR}
onekg_allsamps_dir=${ONEKG_ALLSAMPS_DIR}
onekg_vcf_dir=${ONEKG_VCF_DIR}
onekg_rapid_map_dir=${ONEKG_RAPID_MAP_DIR}
populations=${ONEKG_POPS[*]}
chromosomes=${ONEKG_CHRS[*]}
imaps=${IMAPS[*]}
vcf_tasks=${vcf_tasks}
rapid_map_tasks=${rapid_map_tasks}"

# prepare chromosome-specific population VCFs first
prep_vcf_jid=$(sbatch \
    --parsable \
    --chdir=/home1/karatas/ibd_compare \
    --exclude=d11-05,d11-06,d11-07,b23-11,b23-13,b18-08 \
    --job-name="1kGprepVcf" \
    --array="1-${vcf_tasks}%${MAX_JOBS}" \
    "job_scripts/prepare_1kG_pop_vcfs.sh" \
        "$ONEKG_SOURCE_VCF_DIR" \
        "$ONEKG_ALLSAMPS_DIR" \
        "$ONEKG_VCF_DIR" \
        "${ONEKG_CHRS[@]}"
)
log_msg "submitted 1kG VCF prep; jid=${prep_vcf_jid}"

# prepare rapid maps after successful VCF prep
prep_rapid_map_jid=$(sbatch \
    --parsable \
    --chdir=/home1/karatas/ibd_compare \
    --exclude=d11-05,d11-06,d11-07,b23-11,b23-13,b18-08 \
    --dependency=afterok:${prep_vcf_jid} \
    --job-name="1kGprepMap" \
    --array="1-${rapid_map_tasks}%${MAX_JOBS}" \
    "job_scripts/prepare_1kG_rapid_maps.sh" \
        "$ONEKG_VCF_DIR" \
        "$ONEKG_RAPID_MAP_DIR" \
        "$REF_DIR" \
        "${ONEKG_CHRS[@]}" \
        -- \
        "${IMAPS[@]}" \
        -- \
        "${ONEKG_POPS[@]}"
)
log_msg "submitted 1kG rapid map prep; jid=${prep_rapid_map_jid}; dep=${prep_vcf_jid}"
