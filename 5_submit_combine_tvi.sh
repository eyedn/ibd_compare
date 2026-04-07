#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           5_submit_combine_tvi.sh
###############################################################################


set -euo pipefail
source /home1/karatas/ibd_compare/other_scripts/const.sh
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh

out_agg="${TVI_DIR}/all_stats.tsv"
out_per_rep="${TVI_DIR}/all_stats.per_rep.tsv"

log_msg "Submitting combine_tvi job
tvi_dir=${TVI_DIR}
out_agg=${out_agg}
out_per_rep=${out_per_rep}"

sbatch \
    --job-name="combTVI" \
    --exclude=d11-05,d11-06,d11-07,b23-11,b23-13,b18-08 \
    /home1/karatas/ibd_compare/job_scripts/combine_tvi.sh \
        "$TVI_DIR" \
        "$out_agg" \
        "$out_per_rep"
log_msg "submitted combine_tvi job"
