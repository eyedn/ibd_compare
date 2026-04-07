# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# complete.condition.cols.R
# ______________________________________________________________________________


library(tidyverse)

# standardize condition columns across result tables ----
complete.condition.cols <- function(tbl) {
  required.cols <- c("software.with.params", "imap", "filt", "pop", "threshold")

  for (col.name in required.cols) {
    if (!col.name %in% names(tbl)) {
      tbl[[col.name]] <- NA_character_
    }
  }

  tbl %>%
    select(
      metric,
      question,
      software.with.params,
      imap,
      filt,
      pop,
      threshold,
      contrast.variable,
      contrast,
      level_1,
      level_2,
      level_1.estimate,
      level_2.estimate,
      effect.direction,
      estimate,
      SE,
      df,
      lower.CL,
      upper.CL,
      z.ratio,
      p.value
    ) %>%
    rename(
      contrast.estimate = estimate,
      contrast.se = SE,
      contrast.df = df,
      contrast.ci.lo = lower.CL,
      contrast.ci.hi = upper.CL,
      contrast.z.ratio = z.ratio,
      contrast.p.value = p.value
    )
}
