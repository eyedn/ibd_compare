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

# vcftools command to get files (the --bed is a bed file with intergenic regions)
for f in {1..22}; do 
    vcftools \
        --gzvcf ~/proj2/DataRepository/Human/1000Genome_hg19/vcfs/ALL.chr"$f".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.nodupmarkers.snps.strictmask.vcf.gz \
        --bed chr"$f"_intergenic.bed \
        --keep  ~/proj2/DataRepository/Human/1000GenomeNYGC_hg38/FileInformation/pops/ASW_unrels_no_3rd_N50.txt \
        --het \
        --out chr"$f"
done