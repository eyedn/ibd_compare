# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology 
# Mooney Lab
# ___
# compute.confusion.metrics.windows.R
# ______________________________________________________________________________


library(dplyr)
library(purrr)
library(tibble)
source("job_scripts/R_utils/compute.confusion.metrics.R")

compute.confusion.metrics.windows <- function(
    cov.res,
    inf.ibd,
    region.len.bp,
    n.diploids,
    lower.bounds = c(seq(2, 9.5, by = 0.5), 10),
    upper.bounds = c(seq(2.5, 10, by = 0.5), Inf)
) {
  stopifnot(length(lower.bounds) == length(upper.bounds))

  purrr::map2_dfc(lower.bounds, upper.bounds, function(lo, hi) {

    cov.res.sub <- cov.res %>%
      dplyr::filter(.window.filter(cM.len, lo, hi))

    inf.ibd.sub <- inf.ibd %>%
      dplyr::filter(.window.filter(cM.len, lo, hi))

    # compute metrics on the subset
    out <- compute.confusion.metrics(
      cov.res = cov.res.sub,
      inf.ibd = inf.ibd.sub,
      region.len.bp = region.len.bp,
      n.diploids = n.diploids
    )

    out <- tibble::as_tibble(out)

    # suffix columns with ".<lowerbound>" e.g., tp.2, sensitivity.3.5, ...
    suffix <- format(lo, trim = TRUE, scientific = FALSE)
    names(out) <- paste0(names(out), ".", suffix)

    out
  })
}

# Bins are left-closed except for the final >10 cM bin.
.window.filter <- function(cm.len, lo, hi) {
  if (!is.finite(hi)) {
    return(cm.len > lo)
  }

  if (hi == 10) {
    return(cm.len >= lo & cm.len <= hi)
  }

  cm.len >= lo & cm.len < hi
}
