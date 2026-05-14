# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology 
# Mooney Lab
# ___
# convert_std_bed_to_hapibd.R
# ______________________________________________________________________________


args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(2, 3)) {
  stop(
    "Usage: convert_std_bed_to_hapibd.R <input.bed> <output.gz>",
    " [simulation|empirical]"
  )
}

source("job_scripts/R_utils/std.bed.to.hapibd.R")

fname.in <- args[1]
fname.out <- args[2]
sample.id.mode <- if (length(args) == 3) args[3] else "simulation"

std.bed.to.hapibd(fname.in, fname.out, sample.id.mode)
