#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           const.sh
###############################################################################

# dirs under $HOME
REF_DIR="${HOME}/references"

# dirs under #SCRATCH
OUTDIR="${SCRATCH}/stdpopsim_runs_chr18_admix3" 
COMMON_FF_DIR="${SCRATCH}/ibd_compare_chr18_admix3" 
DEMO_DIR="${SCRATCH}/demo_inf_chr18_admix3"

# dirs under $OUTDIR
TREE_DIR="${OUTDIR}/trees"
VCF_DIR="${OUTDIR}/vcfs"
TRUE_DIR="${OUTDIR}/true_ibd"
ANC_DIR="${OUTDIR}/local_ancestry"
POP_INFO_DIR="${OUTDIR}/pop_info"

# dirs under $COMMON_FF_DIR
TVI_DIR="${COMMON_FF_DIR}/true_vs_inference"

# job array specifications
MAX_JOBS=100

# simulation specifications 
# TMAP exsamples: "HapMapII_GRCh38", "DeCodeSexAveraged_GRCh38", "DefaultUniform"
TMAPS=(
    "HapMapII_GRCh38"
)
POPS=(AFR EUR ADMIX)
NUM_REPS=50
CHR=18
SAMPLE_SIZE=300
MSPRIME_MODEL="dtwf"

# Tennessen 2012 model specifications adpated from Fu 2013
GENERATION_TIME=25
MUTATION_RATE="2.36e-8"
T_AF_YEARS=148000
T_OOA_YEARS=51000
T_EU0_YEARS=23000
T_EG_YEARS=5115
R_EU0="0.00307"
R_EU="0.0195"
R_AF="0.0166"
N_A=7310
N_AF1=14474
N_B=1861
N_EU0=1032
M_AF_B="15e-5"
M_AF_EU="2.5e-5"

# admixture specifications
ADMIXTURE_TIME=15
AFR_PROP=0.75
EUR_PROP=0.25
ADMIX_FOUNDING_SIZE=435
ADMIX_GROWTH_RATE_1="0.02531"
ADMIX_GROWTH_RATE_2="0.01506"
ADMIX_GROWTH_RATE_2_CHANGE_TIME=6
CENSUS_TIME_OFFSET="1e-6"

# ibd ground truth specifications
MIN_SPAN=250000 # measured in bps
MAX_TIME=200 # measured in generations; set to 0 if no max time
MIN_BP_LEN=0
MIN_CM_LEN=2

# ibd inference specifications
IMAPS=(
    "HapMapII_GRCh38"
    "DeCodeSexAveraged_GRCh38"
    "Flat1Mb1cM"
)
SOFTWARE=(
    "hapibd" 
    "rapid"
    )
FILTERING=(
    "none"
    "hot" 
    "ibdends" 
    "hotibdends"
    "ibdendshot"
    )
HOT_BIN_SIZE="0.25"
HOT_BIN_SD=3
KIN_CUTOFF="${KIN_CUTOFF:-0.0442}"

# hapibd params
HAP_MIN_SEED=1.0
HAP_MAX_GAP=2000
HAP_MIN_EXTEND=1.0
HAP_MIN_OUTPUT=2.0
HAP_MIN_MARKERS=100
HAP_MIN_MAC=2

# rapid params
RAP_MIN_LEN=2
RAP_WINDOW_SIZE=300
RAP_NUM_RUNS=10
RAP_NUM_SUCCESSES=2

# demography inference specifications
# gone2 demography inference
FLAT_RATE=1

# ibdne demography inference
