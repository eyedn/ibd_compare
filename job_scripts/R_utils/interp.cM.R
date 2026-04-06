# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology 
# Mooney Lab
# ___
# interp.cM.R
# ______________________________________________________________________________

interp.cM <- function(bp.query, bp.map, cM.map) {
  # limit out-of-bp-range queries
  # below first bp  → 0
  below <- bp.query <= bp.map[1]
  # above last bp   → last cM
  above <- bp.query >= bp.map[length(bp.map)]
  # interior points → interpolate
  mid <- !(below | above)
  
  out <- numeric(length(bp.query))
  out[below] <- 0
  out[above] <- cM.map[length(cM.map)]
  
  # interpolation for mid values
  if (any(mid)) {
    idx.low <- findInterval(bp.query[mid], bp.map)
    
    idx.low[idx.low == 0] <- 1
    idx.low[idx.low == length(bp.map)] <- length(bp.map) - 1
    
    bp.low  <- bp.map[idx.low]
    bp.high <- bp.map[idx.low + 1]
    cM.low  <- cM.map[idx.low]
    cM.high <- cM.map[idx.low + 1]
    
    frac <- (bp.query[mid] - bp.low) / (bp.high - bp.low)
    frac[!is.finite(frac)] <- 0
    
    out[mid] <- cM.low + frac * (cM.high - cM.low)
  }
  
  out
}
  

