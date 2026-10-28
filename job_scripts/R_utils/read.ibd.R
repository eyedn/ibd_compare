# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology 
# Mooney Lab
# ___
# read.ibd.R
# ______________________________________________________________________________


# pattern: Imperative Shell

library(tidyverse)

read.ibd <- function(fname) {
  message("read ibd ", fname)
  # Explicit types also apply when the file contains only its header.
  column.types <- readr::cols(
    chr = readr::col_character(),
    bp.str = readr::col_integer(),
    bp.end = readr::col_integer(),
    bp.len = readr::col_integer(),
    cM.str = readr::col_double(),
    cM.end = readr::col_double(),
    cM.len = readr::col_double(),
    sam1 = readr::col_character(),
    hap1 = readr::col_integer(),
    sam2 = readr::col_character(),
    hap2 = readr::col_integer()
  )
  ibd <- readr::read_tsv(
    fname, col_names = TRUE, col_types = column.types,
    show_col_types = FALSE, lazy = FALSE
  )

  missing.cols <- setdiff(names(column.types$cols), names(ibd))
  if (length(missing.cols) > 0) {
    stop("IBD file ", fname, " is missing standard columns: ",
         paste(missing.cols, collapse = ", "), call. = FALSE)
  }
  if (nrow(readr::problems(ibd)) > 0) {
    stop("Parsing failures in IBD file: ", fname, call. = FALSE)
  }

  ibd
}
