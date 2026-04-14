#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           get_het.sh
###############################################################################

#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=1-00:00:00
#SBATCH --partition=qcb
#SBATCH --account=jazlynmo_738
#SBATCH --nodes=1
#SBATCH --output=/home1/karatas/logs/hapFF.%A_%a.out
#SBATCH --error=/home1/karatas/logs/hapFF.%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=karatas@usc.edu
#SBATCH --array=1-22
#SBATCH --exclude=d11-05,d11-06,d11-07,b23-11,b23-13,b18-08
#SBATCH --chdir=/home1/karatas/ibd_compare
#SBATCH --job-name="ASWhet"

set -euo pipefail
module purge
ml gcc/13.3.0 vcftools/0.1.16
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh

# vcftools command to get files (the --bed is a bed file with intergenic regions)
vcf_dir="${JMPROJ2}/DataRepository/Human/1000GenomeNYGC_hg38/vcfs_strictMask"
intergenic_regions="/home1/karatas/references/GRCh38.p14/gencodeV49.canonicalIntergenic.bed"
ASW_unrels="/home1/karatas/proj/1000GenomeNYGC_RmSamplesGarlic/unrels/ASW.keep"
out_dir="/home1/karatas/proj/1000GenomeNYGC_hg38_karatas/ASW_het"
chrom="${SLURM_ARRAY_TASK_ID}"

mkdir -p "${out_dir}"
log_msg "getting het for chr${chrom}"
vcftools \
    --gzvcf "${vcf_dir}/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chrom}.filtered.shapeit2-duohmm-phased.nodupmarkers.snps.strict.vcf.gz" \
    --bed "$intergenic_regions" \
    --keep "$ASW_unrels" \
    --het \
    --out "${out_dir}/${chrom}"
