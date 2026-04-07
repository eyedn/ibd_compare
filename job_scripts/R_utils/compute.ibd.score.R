# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology 
# Mooney Lab
# ___
# compute.ibd.score.R
# ______________________________________________________________________________


library(tidyverse)

compute.ibd.score <- function(ibd, n.diploids, min.cM = 3, max.cM = 20) {
  message("computing ibd score for this set of ibd, with ", n.diploids, " diploids")
  ibd <- ibd %>% 
    filter(cM.len >= min.cM & cM.len <= max.cM)
  numer <- sum(ibd$bp.len) 
  denom <- choose(2 * n.diploids, 2) - n.diploids
  numer / denom
}