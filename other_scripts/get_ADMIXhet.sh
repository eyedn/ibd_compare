#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           get_ADMIXhet.sh
###############################################################################

#SBATCH --job-name="ADMIXhet"
#SBATCH --chdir=/home1/karatas/ibd_compare
#SBATCH --exclude=d11-05,d11-06,d11-07,b23-11,b23-13,b18-08
#SBATCH --array=1-50
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=1-00:00:00
#SBATCH --partition=qcb
#SBATCH --account=jazlynmo_738
#SBATCH --nodes=1
#SBATCH --output=/home1/karatas/logs/ADMIXhet.%A_%a.out
#SBATCH --error=/home1/karatas/logs/ADMIXhet.%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=karatas@usc.edu


set -euo pipefail
module purge
ml gcc/13.3.0 vcftools/0.1.16
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh

# vcftools command to get files (the --bed is a bed file with intergenic regions)
data_dir="/home1/karatas/scratch/stdpopsim_runs_chr18_admix3"
vcf_dir="${data_dir}/vcfs"
pop_info_dir="${data_dir}/pop_info"
out_dir="${data_dir}/ADMIX_het"

: "${SLURM_ARRAY_TASK_ID:?ERROR: run as a Slurm array}"
vcf="${vcf_dir}/HapMapII_GRCh38_${SLURM_ARRAY_TASK_ID}_ADMIX.vcf.gz"
unrels="${pop_info_dir}/HapMapII_GRCh38_${SLURM_ARRAY_TASK_ID}_ADMIX_unrel_samples.txt"

mkdir -p "${out_dir}"
log_msg "getting ADMIX het for rep${SLURM_ARRAY_TASK_ID}"
vcftools \
    --gzvcf "$vcf" \
    --keep "$unrels" \
    --het \
    --out "${out_dir}/${SLURM_ARRAY_TASK_ID}"
