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
    lower.bounds = c(2, 2.5, 3, 3.5, 4, 4.5, 5),
    upper.bounds = c(2.5, 3, 3.5, 4, 4.5, 5, Inf)
) {
  stopifnot(length(lower.bounds) == length(upper.bounds))

  purrr::map2_dfc(lower.bounds, upper.bounds, function(lo, hi) {

    cov.res.sub <- cov.res %>%
      dplyr::filter(
        cM.len >= lo,
        if (is.finite(hi)) cM.len < hi else TRUE
      )

    inf.ibd.sub <- inf.ibd %>%
      dplyr::filter(
        cM.len >= lo,
        if (is.finite(hi)) cM.len < hi else TRUE
      )

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
