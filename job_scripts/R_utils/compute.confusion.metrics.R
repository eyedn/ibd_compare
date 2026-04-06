# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology 
# Mooney Lab
# ___
# compute.confusion.metrics.R
# ______________________________________________________________________________

library(tidyverse)

compute.confusion.metrics <- function(
    cov.res,
    inf.ibd,
    region.len.bp,
    n.diploids
) {
  
  message("computing confusion matrix metrics for on ", n.diploids, "diploids")
  
  # total possible haplotype pairs in the population
  tot.pos.hap.pairs <- 4 * choose(n.diploids, 2)
  
  # true positive bp: overlapped portion of true IBD segments
  TP <- as.numeric(sum(as.numeric(cov.res$bp.cov), na.rm = TRUE))
  
  # false negative bp: un-covered portion of true IBD segments
  FN <- as.numeric(
    sum(as.numeric(cov.res$bp.len - cov.res$bp.cov), na.rm = TRUE)
    )
  
  # false positive bp: inferred bp that do NOT overlap any true segment
  total.infer.bp <- if (nrow(inf.ibd) == 0) 0 else
    sum(as.numeric(inf.ibd$bp.len), na.rm = TRUE)
  FP <- total.infer.bp - TP
  FP <- as.numeric(pmax(FP, 0))
  
  # true negative: everything else
  total.bp.space <- as.numeric(region.len.bp) * as.numeric(tot.pos.hap.pairs)
  TN <- total.bp.space - TP - FN - FP
  TN <- as.numeric(pmax(TN, 0))
  
  # derive confusion matrix metrics (NA-safe, numeric)
  den_TP_FN <- TP + FN
  den_TN_FP <- TN + FP
  den_TP_FP <- TP + FP
  den_TN_FN <- TN + FN
  den_FP_TN <- FP + TN
  den_ALL   <- TP + TN + FP + FN
  
  sensitivity <- 
    if (!is.na(den_TP_FN) && den_TP_FN > 0) TP / den_TP_FN else NA_real_
  specificity <- 
    if (!is.na(den_TN_FP) && den_TN_FP > 0) TN / den_TN_FP else NA_real_
  PPV <- 
    if (!is.na(den_TP_FP) && den_TP_FP > 0) TP / den_TP_FP else NA_real_
  NPV <- 
    if (!is.na(den_TN_FN) && den_TN_FN > 0) TN / den_TN_FN else NA_real_
  FPR <- 
    if (!is.na(den_FP_TN) && den_FP_TN > 0) FP / den_FP_TN else NA_real_
  accuracy <- 
    if (!is.na(den_ALL)   && den_ALL   > 0) (TP + TN) / den_ALL else NA_real_
  
  data.frame(
    "TP" = as.numeric(TP),
    "FN" = as.numeric(FN),
    "FP" = as.numeric(FP),
    "TN" = as.numeric(TN),
    "sensitivity" = sensitivity,
    "specificity" = specificity,
    "PPV" = PPV,
    "NPV" = NPV,
    "FPR" = FPR,
    "accuracy" = accuracy
    )
}
