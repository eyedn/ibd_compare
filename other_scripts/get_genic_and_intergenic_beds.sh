#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           get_genic_and_intergenic_beds.sh
###############################################################################


set -euo pipefail
module purge
ml gcc/13.3.0 bedtools2/2.31.1

genome_chr_lens="/home1/karatas/references/GRCh38.p14/hg38.chrom.sizes"
canonical_gz="/home1/karatas/references/GRCh38.p14/genesAndGenePredictions.gencodeV49.knownGeneV39.knownCanonical.gz"

gene_bed="/home1/karatas/references/GRCh38.p14/gencodeV49.canonicalGenic.bed"
intergenic_bed="/home1/karatas/references/GRCh38.p14/gencodeV49.canonicalIntergenic.bed"

echo "generate canonical gene bed file; autosomes only"
zcat "$canonical_gz" | \
    awk 'BEGIN{FS="\t"; OFS="\t"}
        nr==1 {next}
        $1 ~ /^chr([1-9]|1[0-9]|2[0-2])$/ {print $1, $2, $3}' | \
    sort -k1,1V -k2,2n | \
    bedtools merge -i - > "$gene_bed"

echo "generate intergenic bed file"
bedtools complement -i "$gene_bed" -g "$genome_chr_lens" > "$intergenic_bed"

echo "done!"
