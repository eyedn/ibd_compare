# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology 
# Mooney Lab
# ___
# count.polytomies.R
# ______________________________________________________________________________

count.polytomies <- function(tr) {
  
  stopifnot(inherits(tr, "phylo"))
  
  n.tips <- length(tr$tip.label)
  n.nodes <- tr$Nnode
  n.total <- n.tips + n.nodes
  
  deg <- integer(n.total)
  for (k in seq_len(nrow(tr$edge))) {
    deg[tr$edge[k, 1]] <- deg[tr$edge[k, 1]] + 1
    deg[tr$edge[k, 2]] <- deg[tr$edge[k, 2]] + 1
  }
  
  internal <- (n.tips + 1):n.total
  
  tibble(
    n.tips = n.tips,
    n.internal = n.nodes,
    n.poly.internal = sum(deg[internal] > 3),
    max.internal.degree = max(deg[internal])
  )
}
