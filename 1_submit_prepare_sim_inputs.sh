#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           1_submit_prepare_sim_inputs.sh
###############################################################################

set -euo pipefail
source /home1/karatas/ibd_compare/other_scripts/const.sh
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh

total_tasks=$(( NUM_REPS * ${#TMAPS[@]} ))
rapid_tasks=$(( total_tasks * ${#IMAPS[@]} * ${#POPS[@]} ))
log_msg "Preparing external OOA simulation inputs: tree_dir=${OOA_TREE_DIR}; vcf_dir=${OOA_VCF_DIR}; outdir=${OUTDIR}; chr=${CHR}; num_reps=${NUM_REPS}; sample_size=${SAMPLE_SIZE}; pops=${POPS[*]}"

prep_inputs_jid=$(sbatch \
    --parsable \
    --chdir=/home1/karatas/ibd_compare \
    --exclude=d11-05,d11-06,d11-07,b23-11,b23-13,b18-08 \
    --job-name="prepIn" \
    --array="1-${total_tasks}%${MAX_JOBS}" \
    job_scripts/prepare_ooa_pop_vcfs.sh \
        "$OOA_VCF_DIR" "$VCF_DIR" "$POP_INFO_DIR" "$CHR" \
        "$SAMPLE_SIZE" "$NUM_REPS" "${TMAPS[@]}" -- "${POPS[@]}"
)
log_msg "submitted OOA VCF preparation; jid=${prep_inputs_jid}"

rapid_jid=$(sbatch \
    --parsable \
    --chdir=/home1/karatas/ibd_compare \
    --exclude=d11-05,d11-06,d11-07,b23-11,b23-13,b18-08 \
    --dependency="afterok:${prep_inputs_jid}" \
    --job-name="prepRAP" \
    --array="1-${rapid_tasks}%${MAX_JOBS}" \
    job_scripts/prepare_rapid_maps.sh \
        "$VCF_DIR" "$RAPID_MAP_DIR" "$REF_DIR" "$CHR" "$NUM_REPS" \
        "${TMAPS[@]}" -- "${IMAPS[@]}" -- "${POPS[@]}"
)
log_msg "submitted RaPID-map preparation; jid=${rapid_jid}"

unrel_jid=$(sbatch \
    --parsable \
    --chdir=/home1/karatas/ibd_compare \
    --exclude=d11-05,d11-06,d11-07,b23-11,b23-13,b18-08 \
    --dependency="afterok:${prep_inputs_jid}" \
    --job-name="unrel" \
    --array="1-${total_tasks}%${MAX_JOBS}" \
    job_scripts/get_unrels.sh "$KIN_CUTOFF" "$VCF_DIR" "$POP_INFO_DIR" \
        "$CHR" "$NUM_REPS" "${TMAPS[@]}" -- "${POPS[@]}"
)
log_msg "submitted unrelated-set extraction; jid=${unrel_jid}"
