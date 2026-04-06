#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           submit_combine_tvi.sh
###############################################################################

set -euo pipefail
source /home1/karatas/ibd_compare/other_scripts/const.sh
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh

OUT_AGG="${TVI_DIR}/all_stats.tsv"
OUT_PER_REP="${TVI_DIR}/all_stats.per_rep.tsv"

log_msg "Submitting combine_tvi job
TVI_DIR=${TVI_DIR}
OUT_AGG=${OUT_AGG}
OUT_PER_REP=${OUT_PER_REP}"

sbatch \
    --job-name="combTVI" \
    --exclude=d11-05,d11-06,d11-07,b23-11,b23-13,b18-08 \
    /home1/karatas/ibd_compare/job_scripts/combine_tvi.sh \
        "$TVI_DIR" \
        "$OUT_AGG" \
        "$OUT_PER_REP"
log_msg "submitted combine_tvi job"
