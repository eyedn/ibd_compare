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

# directory paths
REF_DIR="${HOME}/references"

OUTDIR="${SCRATCH}/stdpopsim_runs_chr18_admix" 
TREE_DIR="${OUTDIR}/trees"
VCF_DIR="${OUTDIR}/vcfs"
TRUE_DIR="${OUTDIR}/true_ibd"
ANC_DIR="${OUTDIR}/local_ancestry"
POP_INFO_DIR="${OUTDIR}/pop_info"

COMMON_FF_DIR="${SCRATCH}/ibd_compare_chr18_admix" 
TVI_DIR="${COMMON_FF_DIR}/true_vs_inference"

DEMO_DIR="${SCRATCH}/demo_inf_chr18_admix"

# job specifications
MAX_JOBS=100

# simulation specifications (e.g., "HapMapII_GRCh38", "DeCodeSexAveraged_GRCh38", "DefaultUniform")
TMAPS=(
    "HapMapII_GRCh38"
)
POPS=(AFR EUR ADMIX)
NUM_REPS=50
# NUM_REPS=2
CHR=18
SAMPLE_SIZE=300 
MSPRIME_MODEL="dtwf"
ADMIXTURE_TIME=14
AFR_PROP=0.75
EUR_PROP=0.25
ADMIX_INITIAL_SIZE=1000
ADMIX_GROWTH_RATE=""
CENSUS_TIME_OFFSET="1e-6"

# ibd specifications
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
