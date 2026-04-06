#!/bin/bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           get_variant_coords.sh
###############################################################################

#SBATCH --job-name=varcoord
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=1-00:00:00
#SBATCH --partition=qcb
#SBATCH --account=jazlynmo_738
#SBATCH --nodes=1
#SBATCH --output=/home1/karatas/logs/slurm.%A_%a.%x.out
#SBATCH --error=/home1/karatas/logs/slurm.%A_%a.%x.err
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=karatas@usc.edu
# NOTE: we'll set the array *at submit time* with --array=...%50
# use... 
#           n=$(ls /home1/karatas/scratch/stdpopsim_runs/*.vcf.gz | wc -l)
#           sbatch --array=0-$((n-1))%100 get_variant_coords.sh

source ~/.bashrc
module purge
module load gcc/13.3.0 htslib/1.19.1 bcftools/1.19

# define paths
VCF_DIR="/home1/karatas/scratch/stdpopsim_runs"
OUT_DIR="${HOME}/scratch/ibd_compare/variant_coords"
mkdir -p "${OUT_DIR}"

# make a bash array of all vcfs
VCF_FILES=("${VCF_DIR}"/*.vcf.gz)

# pick the file for THIS array task
VCF="${VCF_FILES[$SLURM_ARRAY_TASK_ID]}"

# safety: if the index is out of range, just exit
if [[ -z "$VCF" || ! -f "$VCF" ]]; then
    echo "No VCF for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
    exit 0
fi

# build output name: strip .vcf.gz and put in OUT_DIR
base=$(basename "${VCF%.vcf.gz}")
out="${OUT_DIR}/${base}_variant_coords.txt"

bcftools view -H "$VCF" | cut -f1,2 > "$out"
