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
ONEKG_IBD_INIT_INFERENCE_DIR="${SCRATCH}/1kG_ibd_init_inference_ASW_CEU_YRI"
ONEKG_IBD_COMPARE_DIR="${SCRATCH}/1kG_ibd_compare_ASW_CEU_YRI"

# dirs under $OUTDIR
TREE_DIR="${OUTDIR}/trees"
VCF_DIR="${OUTDIR}/vcfs"
PED_DIR="${OUTDIR}/peds"
TRUE_DIR="${OUTDIR}/true_ibd"
ANC_DIR="${OUTDIR}/local_ancestry"
POP_INFO_DIR="${OUTDIR}/pop_info"
RAPID_MAP_DIR="${OUTDIR}/rapid_maps"
ILASH_MAP_DIR="${OUTDIR}/ilash_maps"

# dirs under $COMMON_FF_DIR
TVI_DIR="${COMMON_FF_DIR}/true_vs_inference"

# job array specifications
MAX_JOBS=100
NONADMIXED_CPUS_PER_TASK=4
NONADMIXED_MEM="32G"
NONADMIXED_IBDENDS_THREADS=4
NONADMIXED_IBDENDS_XMX="20G"
ADMIXED_CPUS_PER_TASK=4
ADMIXED_MEM="128G"
ADMIXED_IBDENDS_THREADS=4
ADMIXED_IBDENDS_XMX="96G"

# simulation specifications 
# TMAP exsamples: "HapMapII_GRCh38", "DeCodeSexAveraged_GRCh38", "DefaultUniform"
TMAPS=(
    "HapMapII_GRCh38"
)
POPS=(AFR EUR ADMIX)
NONADMIXED_POP_LABELS=(AFR EUR)
ADMIXED_POP_LABEL="ADMIX"
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
    "ilash"
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
# perturbaton on (HAP_MIN_SEED, HAP_MAX_GAP):
#   (1.0, 2000), (1.0, 1000), (2.0, 1000)
HAP_MIN_SEED=1.0
HAP_MAX_GAP=2000
HAP_MIN_EXTEND=1.0
HAP_MIN_OUTPUT=2.0
HAP_MIN_MARKERS=100
HAP_MIN_MAC=2

# rapid params
# perturbaton on (RAP_WINDOW_SIZE):
#   (100), (200), (300)
RAP_MIN_LEN=2
RAP_WINDOW_SIZE=300
RAP_NUM_RUNS=10
RAP_NUM_SUCCESSES=2

# ilash params
# perturbaton on (ILASH_CM_OVERLAP):
#   (0.5), (1.0), (1.5)
ILASH_AUTO_SLICE=1
ILASH_SLICE_SIZE=350
ILASH_STEP_SIZE=350
ILASH_PERM_COUNT=20
ILASH_SHINGLE_SIZE=15
ILASH_SHINGLE_OVERLAP=0
ILASH_BUCKET_COUNT=5
ILASH_MATCH_THRESHOLD=0.99
ILASH_INTEREST_THRESHOLD=0.70
ILASH_MINHASH_THRESHOLD=55
ILASH_MIN_LENGTH=2.0
ILASH_SLICE_LENGTH=2.0
ILASH_CM_OVERLAP=1.5

# demography inference specifications
# gone2 demography inference
FLAT_RATE=1

# ibdne demography inference

# dirs under $ONEKG_IBD_INIT_INFERENCE_DIR
ONEKG_VCF_DIR="${ONEKG_IBD_INIT_INFERENCE_DIR}/vcfs"
ONEKG_RAPID_MAP_DIR="${ONEKG_IBD_INIT_INFERENCE_DIR}/rapid_maps"

# 1kG data settings
ONEKG_POPS=(ASW CEU YRI)
ONEKG_CHRS=(
    1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 
    )

# read-only 1kG input dirs
ONEKG_SOURCE_VCF_DIR="/home1/karatas/proj/1000GenomeNYGC_hg38_mooney/vcfs_strictMask"
ONEKG_UNREL_DIR="/home1/karatas/proj/1000GenomeNYGC_RmSamplesGarlic/unrels"
ONEKG_ALLSAMPS_DIR="/home1/karatas/proj/1000GenomeNYGC_hg38_mooney/FileInformation/pops/allSamps"
