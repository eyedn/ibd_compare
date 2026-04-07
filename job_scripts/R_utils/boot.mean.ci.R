# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology 
# Mooney Lab
# ___
# boot.mean.ci.R
# ______________________________________________________________________________


boot.mean.ci <- function(x, R = 1000, conf = 0.95) {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    return(c(mean = NA_real_, lo = NA_real_, hi = NA_real_))
  }
  n <- length(x)
  boots <- replicate(R, mean(sample(x, n, replace = TRUE)))
  alpha <- (1 - conf) / 2
  c(
    mean = mean(x),
    lo = unname(quantile(boots, alpha)),
    hi = unname(quantile(boots, 1 - alpha))
  )
}