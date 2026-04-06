# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology 
# Mooney Lab
# ___
# read.cov.res.R
# ______________________________________________________________________________

library(tidyverse)

read.cov.res <- function(fname) {
  message("read coverage results from  ", fname)
  
  cov.res <- readr::read_tsv(
    file = fname, col_names = FALSE, show_col_types = FALSE
    )
  
  # preferred specific formatting
  cov.res %>%
    rename(
      chr = X1, 
      bp.str = X2, 
      bp.end = X3, 
      bp.len = X4, 
      cM.str = X5, 
      cM.end = X6, 
      cM.len = X7, 
      sam1 = X8, 
      hap1 = X9, 
      sam2 = X10, 
      hap2 = X11, 
      num.infer.ibd = X12, 
      bp.cov = X13, 
      bp.perc.cov = X14
      ) %>%
    mutate(
      chr = as.character(chr), 
      bp.str = as.integer(bp.str), 
      bp.end = as.integer(bp.end), 
      bp.len = as.integer(bp.len), 
      cM.str = as.numeric(cM.str), 
      cM.end = as.numeric(cM.end), 
      cM.len = as.numeric(cM.len), 
      sam1 = as.integer(sam1), 
      hap1 = as.integer(hap1), 
      sam2 = as.integer(sam2), 
      hap2 = as.integer(hap2), 
      num.infer.ibd = as.integer(num.infer.ibd), 
      bp.cov = as.integer(bp.cov), 
      bp.perc.cov = as.numeric(bp.perc.cov)
    )
  
}
