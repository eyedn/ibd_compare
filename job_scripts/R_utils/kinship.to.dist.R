# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology 
# Mooney Lab
# ___
# kinship.to.dist.R
# ______________________________________________________________________________

library(tidyverse)

kinship.to.dist <- function(kin.mat, eps = 0) {
  
  message("converting kinship to distance...")
  
  stopifnot(is.matrix(kin.mat), nrow(kin.mat) == ncol(kin.mat))
  
  kin.mat <- (kin.mat + t(kin.mat)) / 2
  
  # make sure diagonal is defined (convention: 0.5)
  if (anyNA(diag(kin.mat))) diag(kin.mat) <- 0.5
  
  # d_ij = sqrt(Kii + Kjj - 2Kij)
  kii  <- as.numeric(diag(kin.mat))
  d.mat <- outer(kii, kii, "+") - 2 * kin.mat
  d.mat <- pmax(d.mat, 0)
  d.mat <- sqrt(d.mat)          
  diag(d.mat) <- 0
  
  # optional: tiny jitter to break ties deterministically
  if (!is.null(eps) && eps > 0) {
    set.seed(1)
    j <- matrix(runif(length(d.mat), -eps, eps), nrow = nrow(d.mat))
    d.mat <- d.mat + (j + t(j)) / 2
    diag(d.mat) <- 0
    d.mat <- pmax(d.mat, 0)
  }
  
  as.dist(d.mat)
}
