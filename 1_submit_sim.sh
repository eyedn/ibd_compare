#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           1_submit_sim.sh
###############################################################################


set -euo pipefail
source /home1/karatas/ibd_compare/other_scripts/const.sh
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh

# NOTE:  sim. specs in /home1/karatas/ibd_compare/const.sh
total_tasks=$(( NUM_REPS * ${#TMAPS[@]} ))
log_msg "Detected sim specifications
outdir=${OUTDIR}
num_reps=${NUM_REPS}
sample_size=${SAMPLE_SIZE}
tree_dir=${TREE_DIR}
msprime_model=${MSPRIME_MODEL}
anc_dir=${ANC_DIR}
chr=${CHR}
total_tasks=${total_tasks} (per job array)"

# submit sims
total_tasks=$(( NUM_REPS * ${#TMAPS[@]} ))
sim_jid=$(sbatch \
    --parsable \
    --chdir=/home1/karatas/ibd_compare \
    --exclude=d11-05,d11-06,d11-07,b23-11,b23-13,b18-08 \
    --job-name="sim" \
    --array="1-${total_tasks}%${MAX_JOBS}" \
    "job_scripts/sim_ooaa.sh" \
        "$TREE_DIR" "$SAMPLE_SIZE" "$NUM_REPS" \
        "$MSPRIME_MODEL" "$CHR" \
        "$GENERATION_TIME" "$MUTATION_RATE" \
        "$T_AF_YEARS" "$T_OOA_YEARS" "$T_EU0_YEARS" "$T_EG_YEARS" \
        "$R_EU0" "$R_EU" "$R_AF" \
        "$N_A" "$N_AF1" "$N_B" "$N_EU0" \
        "$M_AF_B" "$M_AF_EU" \
        "$ADMIXTURE_TIME" "$AFR_PROP" "$EUR_PROP" \
        "$ADMIX_INITIAL_SIZE" "$ADMIX_GROWTH_RATE_1" \
        "$ADMIX_GROWTH_RATE_2" "$ADMIX_GROWTH_RATE_2_CHANGE_TIME" \
        "$CENSUS_TIME_OFFSET" \
        "${TMAPS[@]}"
)
log_msg "submitted combined sim array; jid=${sim_jid}"

anc_jid=$(sbatch \
    --parsable \
    --chdir=/home1/karatas/ibd_compare \
    --exclude=d11-05,d11-06,d11-07,b23-11,b23-13,b18-08 \
    --dependency=afterok:${sim_jid} \
    --job-name="anc" \
    --array="1-${total_tasks}%${MAX_JOBS}" \
    "job_scripts/get_local_ancestry.sh" \
        "$TREE_DIR" \
        "$ANC_DIR" \
        "$CHR" \
        "$SAMPLE_SIZE" \
        "$ADMIXTURE_TIME" \
        "$CENSUS_TIME_OFFSET" \
        "$NUM_REPS" \
        "${TMAPS[@]}" \
        -- \
        "${POPS[@]}"
)
log_msg "local ancestry array for each pop and map; jid=${anc_jid}; dep on ${sim_jid}"

# subset vcf file by population
total_tasks=$(( NUM_REPS * ${#TMAPS[@]} ))
subset_jid=$(sbatch \
    --parsable \
    --chdir=/home1/karatas/ibd_compare \
    --exclude=d11-05,d11-06,d11-07,b23-11,b23-13,b18-08 \
    --dependency=afterok:${sim_jid} \
    --job-name="sub" \
    --array="1-${total_tasks}%${MAX_JOBS}" \
    "job_scripts/subset_ts.sh" \
        "$TREE_DIR" \
        "$VCF_DIR" \
        "$POP_INFO_DIR" \
        "$CHR" \
        "$SAMPLE_SIZE" \
        "$NUM_REPS" \
        "${TMAPS[@]}" \
        -- \
        "${POPS[@]}"
)
log_msg "subset ts array for each pop and create vcfs; jid=${subset_jid}; dep on ${sim_jid}"

# get unrelated individuals from subsetted vcf files
unrel_jid=$(sbatch \
    --parsable \
    --chdir=/home1/karatas/ibd_compare \
    --exclude=d11-05,d11-06,d11-07,b23-11,b23-13,b18-08 \
    --dependency=afterok:${subset_jid} \
    --job-name="unrel" \
    --array="1-${total_tasks}%${MAX_JOBS}" \
    "job_scripts/get_unrels.sh" \
        "$KIN_CUTOFF" \
        "$VCF_DIR" \
        "$POP_INFO_DIR" \
        "$NUM_REPS" \
        "${TMAPS[@]}" \
        -- \
        "${POPS[@]}"
)
log_msg "subset vcf unrels for each pop and map; jid=${unrel_jid}; dep. on ${subset_jid}"
