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
  tp <- as.numeric(sum(as.numeric(cov.res$bp.cov), na.rm = TRUE))
  
  # false negative bp: un-covered portion of true IBD segments
  fn <- as.numeric(
    sum(as.numeric(cov.res$bp.len - cov.res$bp.cov), na.rm = TRUE)
    )
  
  # false positive bp: inferred bp that do NOT overlap any true segment
  total.infer.bp <- if (nrow(inf.ibd) == 0) 0 else
    sum(as.numeric(inf.ibd$bp.len), na.rm = TRUE)
  fp <- total.infer.bp - tp
  fp <- as.numeric(pmax(fp, 0))
  
  # true negative: everything else
  total.bp.space <- as.numeric(region.len.bp) * as.numeric(tot.pos.hap.pairs)
  tn <- total.bp.space - tp - fn - fp
  tn <- as.numeric(pmax(tn, 0))
  
  # derive confusion matrix metrics (NA-safe, numeric)
  den.tp.fn <- tp + fn
  den.tn.fp <- tn + fp
  den.tp.fp <- tp + fp
  den.tn.fn <- tn + fn
  den.fp.tn <- fp + tn
  den.all   <- tp + tn + fp + fn
  
  sensitivity <- 
    if (!is.na(den.tp.fn) && den.tp.fn > 0) tp / den.tp.fn else NA_real_
  specificity <- 
    if (!is.na(den.tn.fp) && den.tn.fp > 0) tn / den.tn.fp else NA_real_
  ppv <- 
    if (!is.na(den.tp.fp) && den.tp.fp > 0) tp / den.tp.fp else NA_real_
  npv <- 
    if (!is.na(den.tn.fn) && den.tn.fn > 0) tn / den.tn.fn else NA_real_
  fpr <- 
    if (!is.na(den.fp.tn) && den.fp.tn > 0) fp / den.fp.tn else NA_real_
  accuracy <- 
    if (!is.na(den.all)   && den.all   > 0) (tp + tn) / den.all else NA_real_
  
  data.frame(
    "tp" = as.numeric(tp),
    "fn" = as.numeric(fn),
    "fp" = as.numeric(fp),
    "tn" = as.numeric(tn),
    "sensitivity" = sensitivity,
    "specificity" = specificity,
    "ppv" = ppv,
    "npv" = npv,
    "fpr" = fpr,
    "accuracy" = accuracy
    )
}
