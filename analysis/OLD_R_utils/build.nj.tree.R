# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology 
# Mooney Lab
# ___
# build.nj.tree.R
# ______________________________________________________________________________

library(ape)

# dist -> NJ tree, then collapse short edges to allow polytomies
build.nj.tree <- function(kin.mat, collapse.len = 0, jitter.eps = 0) {
  
  message("building NJ tree...")
  
  d <- kinship.to.dist(kin.mat, eps = jitter.eps)
  tr <- ape::nj(d)
  tr <- ape::ladderize(tr)
  
  if (!is.null(collapse.len) && collapse.len > 0) {
    message("collapsing edges <= ", collapse.len, " to allow polytomies...")
    tr$edge.length[is.na(tr$edge.length)] <- 0
    tr <- ape::di2multi(tr, tol = collapse.len)
  }
  
  tr
}