# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# build.confusion.model.formula.R
# ______________________________________________________________________________


# build a confusion-metric model formula from included fixed terms ----
build.confusion.model.formula <- function(fixed.terms) {
  stopifnot(length(fixed.terms) > 0)

  stats::as.formula(
    paste(
      "cbind(successes, failures) ~",
      paste(fixed.terms, collapse = " + "),
      "+ (1 | replicate.id)"
    )
  )
}
