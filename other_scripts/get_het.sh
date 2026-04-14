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


set -euo pipefail
module purge
ml gcc/13.3.0 vcftools/0.1.16
source /home1/karatas/ibd_compare/other_scripts/log_msg.sh

# vcftools command to get files (the --bed is a bed file with intergenic regions)
vcf_dir="${JMPROJ2}/DataRepository/Human/1000GenomeNYGC_hg38/vcfs_strictMask"
intergenic_regions="/home1/karatas/references/GRCh38.p14/gencodeV49.canonicalIntergenic.bed"
ASW_unrels="/home1/karatas/proj/1000GenomeNYGC_RmSamplesGarlic/unrels/ASW.keep"
out_dir="/home1/karatas/proj/1000GenomeNYGC_hg38_karatas/ASW_het"

mkdir -p "${out_dir}"
for f in {1..22}; do 
    log_msg "getting het for chr${f}"
    vcftools \
        --gzvcf "${vcf_dir}/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${f}.filtered.shapeit2-duohmm-phased.nodupmarkers.snps.strict.vcf.gz" \
        --bed "$intergenic_regions" \
        --keep "$ASW_unrels" \
        --het \
        --out "${out_dir}/${f}"
done