#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           get_variant_coords.sh
###############################################################################


source ~/.bashrc
module purge
module load gcc/13.3.0 htslib/1.19.1 bcftools/1.19

# define paths
vcf_dir="/home1/karatas/scratch/stdpopsim_runs"
out_dir="${HOME}/scratch/ibd_compare/variant_coords"
mkdir -p "${out_dir}"

# make a bash array of all vcfs
vcf_files=("${vcf_dir}"/*.vcf.gz)

# pick the file for THIS array task
vcf="${vcf_files[$SLURM_ARRAY_TASK_ID]}"

# safety: if the index is out of range, just exit
if [[ -z "$vcf" || ! -f "$vcf" ]]; then
    echo "No vcf for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
    exit 0
fi

# build output name: strip .vcf.gz and put in out_dir
base=$(basename "${vcf%.vcf.gz}")
out="${out_dir}/${base}_variant_coords.txt"

bcftools view -H "$vcf" | cut -f1,2 > "$out"
