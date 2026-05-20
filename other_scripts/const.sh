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

# simulation specifications: 4B18 or 3T
MODEL_SUFFIX="${MODEL_SUFFIX:-3T}"

# TMAP exsamples: "HapMapII_GRCh38", "DeCodeSexAveraged_GRCh38", "DefaultUniform"
TMAPS=(
    "HapMapII_GRCh38"
)
NUM_REPS=50
CHR=18
SAMPLE_SIZE=300
MSPRIME_MODEL="dtwf"
STDPOPSIM_DEMOGRAPHY_ID="AmericanAdmixture_4B18"

case "${MODEL_SUFFIX}" in
    3T)
        POPS=(AFR EUR ADMIX)
        NONADMIXED_POP_LABELS=(AFR EUR)
        ADMIXED_POP_LABEL="ADMIX"
        ;;
    4B18)
        POPS=(AFR EUR ASIA ADMIX)
        NONADMIXED_POP_LABELS=(AFR EUR ASIA)
        ADMIXED_POP_LABEL="ADMIX"
        ;;
    *)
        echo "ERROR: MODEL_SUFFIX must be 3T or 4B18" >&2
        return 1 2>/dev/null || exit 1
        ;;
esac

# dirs under $HOME
REF_DIR="${HOME}/references"

# dirs under #SCRATCH
OUTDIR="${SCRATCH}/stdpopsim_runs_chr${CHR}_admix_${MODEL_SUFFIX}"
COMMON_FF_DIR="${SCRATCH}/ibd_compare_chr${CHR}_admix_${MODEL_SUFFIX}"
DEMO_DIR="${SCRATCH}/demo_inf_chr${CHR}_admix_${MODEL_SUFFIX}"
ONEKG_IBD_INIT_INFERENCE_DIR="${SCRATCH}/1kG_ibd_init_inference_1kG"
ONEKG_IBD_COMPARE_DIR="${SCRATCH}/1kG_ibd_compare_1kG"

# dirs under $OUTDIR
TREE_DIR="${OUTDIR}/trees"
VCF_DIR="${OUTDIR}/vcfs"
PED_DIR="${OUTDIR}/peds"
TRUE_DIR="${OUTDIR}/true_ibd"
ANC_DIR="${OUTDIR}/local_ancestry"
GLOBAL_ANC_DIR="${OUTDIR}/global_ancestry"
POP_INFO_DIR="${OUTDIR}/pop_info"
RAPID_MAP_DIR="${OUTDIR}/rapid_maps"
ILASH_MAP_DIR="${OUTDIR}/ilash_maps"

# dirs under $COMMON_FF_DIR
TVI_DIR="${COMMON_FF_DIR}/true_vs_inference"

# job array specifications
MAX_JOBS=100
NONADMIXED_CPUS_PER_TASK=1
NONADMIXED_MEM="128G"
NONADMIXED_IBDENDS_THREADS=1
NONADMIXED_IBDENDS_XMX="96G"
ADMIXED_CPUS_PER_TASK=4
ADMIXED_MEM="128G"
ADMIXED_IBDENDS_THREADS=4
ADMIXED_IBDENDS_XMX="96G"

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
CENSUS_TIME_OFFSET="1e-6"
ADMIX_GENERATION_COUNT=16
ADMIX_MIXING_GENERATION_COUNT=11
ADMIX_FOUNDER_AFR_COUNT=183
ADMIX_FOUNDER_EUR_COUNT=25
ADMIX_MODERN_GROWTH_RATE="0.023175"
ADMIX_NE_BY_GENERATION=(
    187.6667
    2.214286
    276.3148
    2019.588
    10015.88
    45816.03
    149343.0
    381940.1
    390895.0
    400059.7
    409439.4
    419039.0
    428863.6
    438918.6
    449209.3
    459741.3
)
ADMIX_AFR_PROPS_BY_GENERATION=(
    0.8798077
    0.920000
    0.904820
    0.791384
    0.786692
    0.719992
    0.494960
    0.130456
    0.060000
    0.060000
    0.060000
    0.000000
    0.000000
    0.000000
    0.000000
    0.000000
)
ADMIX_EUR_PROPS_BY_GENERATION=(
    0.1201923
    0.080000
    0.080000
    0.080000
    0.080000
    0.080000
    0.080000
    0.080000
    0.030000
    0.030000
    0.030000
    0.000000
    0.000000
    0.000000
    0.000000
    0.000000
)
ADMIX_HYBRID_PROPS_BY_GENERATION=(
    0.000000
    0.000000
    0.015180
    0.128616
    0.133308
    0.200008
    0.425040
    0.789544
    0.910000
    0.910000
    0.910000
    1.000000
    1.000000
    1.000000
    1.000000
    1.000000
)

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

# ibd-ends params
IBDENDS_MIN_MAF=0.001

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
