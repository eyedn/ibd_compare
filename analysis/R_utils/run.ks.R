# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# run.ks.R
# ______________________________________________________________________________


run.ks <- function(tru, inf) {
  tru <- tru[is.finite(tru)]
  inf <- inf[is.finite(inf)]
  if (length(tru) < 10 || length(inf) < 10) return(NA_real_)
  ks.test(tru, inf)$p.value
}
