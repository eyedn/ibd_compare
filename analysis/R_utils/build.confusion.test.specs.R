# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# build.confusion.test.specs.R
# ______________________________________________________________________________


library(tidyverse)

# define the omnibus lrt terms for overall and bin-effect models ----
build.confusion.test.specs <- function(model.scope) {
  stopifnot(model.scope %in% c("overall", "bin"))

  if (model.scope == "overall") {
    full.terms <- c("software.with.params", "imap", "filt", "pop")
    tested.terms <- full.terms
  } else {
    full.terms <- c(
      "software.with.params",
      "imap",
      "filt",
      "pop",
      "threshold",
      "software.with.params:threshold",
      "imap:threshold",
      "filt:threshold",
      "pop:threshold"
    )
    tested.terms <- c(
      "threshold",
      "software.with.params:threshold",
      "imap:threshold",
      "filt:threshold",
      "pop:threshold"
    )
  }

  tibble(
    model.scope = model.scope,
    tested.term = tested.terms,
    reduced.terms = purrr::map(
      tested.terms,
      ~ setdiff(full.terms, .x)
    ),
    full.terms = list(full.terms)
  )
}
