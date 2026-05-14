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

std.bed.to.hapibd <- function(fname.in, fname.out,
                              sample.id.mode = c("simulation", "empirical")) {
  sample.id.mode <- match.arg(sample.id.mode)
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

  # Keep the simulation-only tsk_* contract separate from empirical labels.
  format.sample.ids <- function(sample_ids) {
    sample_ids <- as.character(sample_ids)

    if (sample.id.mode == "empirical") {
      return(sample_ids)
    }

    sample_ids_int <- suppressWarnings(as.integer(sample_ids))
    if (any(is.na(sample_ids_int))) {
      stop(
        "simulation sample.id.mode requires integer-compatible sample IDs"
      )
    }

    paste0("tsk_", sample_ids_int)
  }

  hapibd.calls <- calls %>%
    arrange(bp.str, bp.end) %>%
    transmute(
      sam1 = format.sample.ids(sam1),
      hap1 = as.integer(hap1),
      sam2 = format.sample.ids(sam2),
      hap2 = as.integer(hap2),
      chr = chr,
      bp.str = as.integer(bp.str),
      bp.end = as.integer(bp.end),
      cM.len = as.numeric(cM.len)
    )

  readr::write_tsv(hapibd.calls, fname.out, col_names = FALSE)
}
