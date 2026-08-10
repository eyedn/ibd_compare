# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology 
# Mooney Lab
# ___
# std.bed.to.hapibd.R
# ______________________________________________________________________________


library(tidyverse)

std.bed.to.hapibd <- function(fname.in, fname.out) {
  req.cols <- c(
    "chr", "bp.str", "bp.end", "bp.len", "cM.str", "cM.end", "cM.len",
    "sam1", "hap1", "sam2", "hap2"
  )

  calls <- readr::read_tsv(
    fname.in, col_names = TRUE, show_col_types = FALSE
  )

  missing.cols <- setdiff(req.cols, colnames(calls))
  if (length(missing.cols) > 0) {
    stop(
      "std BED file missing required columns: ",
      paste(missing.cols, collapse = ", ")
    )
  }

  hapibd.calls <- calls %>%
    arrange(bp.str, bp.end) %>%
    transmute(
      sam1 = as.character(sam1),
      hap1 = as.integer(hap1),
      sam2 = as.character(sam2),
      hap2 = as.integer(hap2),
      chr = chr,
      bp.str = as.integer(bp.str),
      bp.end = as.integer(bp.end),
      cM.len = as.numeric(cM.len)
    )

  readr::write_tsv(hapibd.calls, fname.out, col_names = FALSE)
}
