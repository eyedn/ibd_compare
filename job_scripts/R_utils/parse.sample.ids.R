# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# parse.sample.ids.R
# ______________________________________________________________________________


library(tidyverse)

parse.sample.ids <- function(sample_ids, sample.id.mode,
                             non.integer.error.message) {
  sample_ids <- as.character(sample_ids)

  if (sample.id.mode == "empirical") {
    return(sample_ids)
  }

  sample_ids <- str_extract(sample_ids, "\\d+$")
  sample_ids <- suppressWarnings(as.integer(sample_ids))

  if (any(is.na(sample_ids))) {
    stop(non.integer.error.message)
  }

  sample_ids
}
