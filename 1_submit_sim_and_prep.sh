#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           1_submit_sim_and_prep.sh
###############################################################################


set -euo pipefail
source /home1/karatas/ibd_compare/other_scripts/const.sh
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh

# NOTE:  sim. specs in /home1/karatas/ibd_compare/const.sh
total_tasks=$(( NUM_REPS * ${#TMAPS[@]} ))
prep_total_tasks=$(( NUM_REPS * ${#TMAPS[@]} * ${#IMAPS[@]} * ${#POPS[@]} ))
log_msg "Detected sim+prep specifications
model_suffix=${MODEL_SUFFIX}
outdir=${OUTDIR}
num_reps=${NUM_REPS}
sample_size=${SAMPLE_SIZE}
tree_dir=${TREE_DIR}
msprime_model=${MSPRIME_MODEL}
anc_dir=${ANC_DIR}
ped_dir=${PED_DIR}
rapid_map_dir=${RAPID_MAP_DIR}
ilash_map_dir=${ILASH_MAP_DIR}
chr=${CHR}
total_tasks=${total_tasks} (per job array)
prep_map_tasks=${prep_total_tasks} (combined rapid+ilash map tasks)
pops=${POPS[*]}"

# submit sims
case "${MODEL_SUFFIX}" in
    3T)
        sim_jid=$(sbatch \
            --parsable \
            --chdir=/home1/karatas/ibd_compare \
            --exclude=d11-05,d11-06,d11-07,b23-11,b23-13,b18-08 \
            --job-name="sim" \
            --array="1-${total_tasks}%${MAX_JOBS}" \
            "job_scripts/sim_3T.sh" \
                "$TREE_DIR" "$SAMPLE_SIZE" "$NUM_REPS" \
                "$MSPRIME_MODEL" "$CHR" \
                "$GENERATION_TIME" "$MUTATION_RATE" \
                "$T_AF_YEARS" "$T_OOA_YEARS" "$T_EU0_YEARS" "$T_EG_YEARS" \
                "$R_EU0" "$R_EU" "$R_AF" \
                "$N_A" "$N_AF1" "$N_B" "$N_EU0" \
                "$M_AF_B" "$M_AF_EU" \
                "$ADMIXTURE_TIME" "$AFR_PROP" "$EUR_PROP" \
                "$ADMIX_FOUNDING_SIZE" "$ADMIX_GROWTH_RATE_1" \
                "$ADMIX_GROWTH_RATE_2" "$ADMIX_GROWTH_RATE_2_CHANGE_TIME" \
                "$CENSUS_TIME_OFFSET" \
                "${TMAPS[@]}"
        )
        log_msg "submitted 3T sim array via sim_3T.sh; jid=${sim_jid}"

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
        log_msg "submitted local ancestry array; jid=${anc_jid}; dep on ${sim_jid}"
        ;;
    4B18)
        sim_jid=$(sbatch \
            --parsable \
            --chdir=/home1/karatas/ibd_compare \
            --exclude=d11-05,d11-06,d11-07,b23-11,b23-13,b18-08 \
            --job-name="sim" \
            --array="1-${total_tasks}%${MAX_JOBS}" \
            "job_scripts/sim_4B18.sh" \
                "$TREE_DIR" \
                "$SAMPLE_SIZE" \
                "$NUM_REPS" \
                "$MSPRIME_MODEL" \
                "$STDPOPSIM_DEMOGRAPHY_ID" \
                "$CHR" \
                "${TMAPS[@]}"
        )
        log_msg "submitted 4B18 sim array via sim_4B18.sh; jid=${sim_jid}"
        log_msg "skipping local ancestry for MODEL_SUFFIX=4B18"
        ;;
esac

# prepare population-wise vcfs and peds
prep_inputs_jid=$(sbatch \
    --parsable \
    --chdir=/home1/karatas/ibd_compare \
    --exclude=d11-05,d11-06,d11-07,b23-11,b23-13,b18-08 \
    --dependency=afterok:${sim_jid} \
    --job-name="prepIn" \
    --array="1-${total_tasks}%${MAX_JOBS}" \
    "job_scripts/prepare_pop_vcfs_peds.sh" \
        "$TREE_DIR" \
        "$VCF_DIR" \
        "$PED_DIR" \
        "$POP_INFO_DIR" \
        "$CHR" \
        "$SAMPLE_SIZE" \
        "$NUM_REPS" \
        "${TMAPS[@]}" \
        -- \
        "${POPS[@]}"
)
log_msg "prep vcfs+peds array for each pop and map; jid=${prep_inputs_jid}; dep on ${sim_jid}"

# prepare rapid+ilash maps for each tmap x imap x seed x pop
prep_map_jid=$(sbatch \
    --parsable \
    --chdir=/home1/karatas/ibd_compare \
    --exclude=d11-05,d11-06,d11-07,b23-11,b23-13,b18-08 \
    --dependency=afterok:${prep_inputs_jid} \
    --job-name="prepRIL" \
    --array="1-${prep_total_tasks}%${MAX_JOBS}" \
    "job_scripts/prepare_rapid_ilash_maps.sh" \
        "$VCF_DIR" \
        "$PED_DIR" \
        "$RAPID_MAP_DIR" \
        "$ILASH_MAP_DIR" \
        "$REF_DIR" \
        "$CHR" \
        "$NUM_REPS" \
        "${TMAPS[@]}" \
        -- \
        "${IMAPS[@]}" \
        -- \
        "${POPS[@]}"
)
log_msg "prep rapid+ilash maps array; jid=${prep_map_jid}; dep on ${prep_inputs_jid}"

# get unrelated individuals from subsetted vcf files
unrel_jid=$(sbatch \
    --parsable \
    --chdir=/home1/karatas/ibd_compare \
    --exclude=d11-05,d11-06,d11-07,b23-11,b23-13,b18-08 \
    --dependency=afterok:${prep_inputs_jid} \
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
log_msg "subset vcf unrels for each pop and map; jid=${unrel_jid}; dep. on ${prep_inputs_jid}"
