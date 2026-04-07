# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology 
# Mooney Lab
# ___
# compute.gen.rf.R
# ______________________________________________________________________________


library(ape)
library(tidyverse)
library(TreeDist)

compute.gen.rf <- function(tree.true, tree.inf) {
  
  message("computing generalized RF distance (TreeDist)...")
  
  stopifnot(
    inherits(tree.true, "phylo"),
    inherits(tree.inf, "phylo"),
    setequal(tree.true$tip.label, tree.inf$tip.label)
  )
  
  # enforce identical tip order
  tree.inf <- ape::keep.tip(tree.inf, tree.true$tip.label)
  
  # recommended by TreeDist docs (handles polytomies sensibly)
  d <- TreeDist::ClusteringInfoDistance(tree.true, tree.inf)
  
  tibble(gen.rf = as.numeric(d))
}
