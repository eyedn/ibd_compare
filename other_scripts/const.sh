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

# Read-only OOA simulation inputs.
OOA_TREE_DIR="${SCRATCH}/OOA_NAAdmixture_small/trees"
OOA_VCF_DIR="${SCRATCH}/OOA_NAAdmixture_small/vcfs"
TMAPS=("HapMapII_GRCh38")
POPS=(AFR EUR ADX)
NONADMIXED_POP_LABELS=(AFR EUR)
ADMIXED_POP_LABEL="ADX"
NUM_REPS=50
SAMPLE_SIZE=500
CHR=18

# IBD-comparison outputs are deliberately separate from OOA inputs.
REF_DIR="${HOME}/references"
OUTDIR="${SCRATCH}/ibd_compare_chr${CHR}_ooa"
VCF_DIR="${OUTDIR}/vcfs"
POP_INFO_DIR="${OUTDIR}/pop_info"
RAPID_MAP_DIR="${OUTDIR}/rapid_maps"
TRUE_DIR="${OUTDIR}/true_ibd"
COMMON_FF_DIR="${OUTDIR}/filter_format"
TVI_DIR="${COMMON_FF_DIR}/true_vs_inference"
ONEKG_IBD_INIT_INFERENCE_DIR="${SCRATCH}/ibd_compare_init_1kG"
ONEKG_IBD_COMPARE_DIR="${SCRATCH}/ibd_compare_1kG"

# Job-array and filtering resources.
MAX_JOBS=100
NONADMIXED_CPUS_PER_TASK=4
NONADMIXED_MEM="128G"
NONADMIXED_IBDENDS_THREADS=4
NONADMIXED_IBDENDS_XMX="96G"
ADMIXED_CPUS_PER_TASK=4
ADMIXED_MEM="128G"
ADMIXED_IBDENDS_THREADS=4
ADMIXED_IBDENDS_XMX="96G"

# IBD ground-truth specifications.
MAX_TIME=75
MIN_SPAN=100000
MIN_BP_LEN=0
MIN_CM_LEN=2

# IBD inference specifications.
IMAPS=("HapMapII_GRCh38" "DeCodeSexAveraged_GRCh38" "Flat1Mb1cM")
SOFTWARE=("hapibd" "rapid")
FILTERING=("none" "hot" "ibdends" "hotibdends" "ibdendshot")
HOT_BIN_SIZE="0.25"
HOT_BIN_SD=3
KIN_CUTOFF="${KIN_CUTOFF:-0.0442}"

# hap-IBD parameters.
HAP_MIN_SEED=2.0
HAP_MAX_GAP=1000
HAP_MIN_EXTEND=1.0
HAP_MIN_OUTPUT=2.0
HAP_MIN_MARKERS=100
HAP_MIN_MAC=2

# RaPID parameters.
RAP_MIN_LEN=2
RAP_WINDOW_SIZE=100
RAP_NUM_RUNS=10
RAP_NUM_SUCCESSES=2

# ibd-ends parameters.
IBDENDS_MIN_MAF=0.001

# 1kG output directories and read-only inputs.
ONEKG_VCF_DIR="${ONEKG_IBD_INIT_INFERENCE_DIR}/vcfs"
ONEKG_RAPID_MAP_DIR="${ONEKG_IBD_INIT_INFERENCE_DIR}/rapid_maps"
ONEKG_POPS=(ASW CEU YRI)
ONEKG_CHRS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22)
ONEKG_SOURCE_VCF_DIR="/home1/karatas/proj/1000GenomeNYGC_hg38_mooney/vcfs_strictMask"
ONEKG_UNREL_DIR="/home1/karatas/proj/1000GenomeNYGC_RmSamplesGarlic/unrels"
ONEKG_ALLSAMPS_DIR="/home1/karatas/proj/1000GenomeNYGC_hg38_mooney/FileInformation/pops/allSamps"
