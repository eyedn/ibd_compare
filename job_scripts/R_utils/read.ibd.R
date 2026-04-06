# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology 
# Mooney Lab
# ___
# read.ibd.R
# ______________________________________________________________________________

library(tidyverse)

read.ibd <- function(fname) {
  message("read ibd ", fname)
  readr::read_tsv(fname, col_names = TRUE, show_col_types = FALSE)
}