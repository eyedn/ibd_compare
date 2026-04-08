# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# tidy.confusion.emmeans.R
# ______________________________________________________________________________


library(tidyverse)
library(emmeans)

# extract response-scale marginal means or pairwise contrasts ----
tidy.confusion.emmeans <- function(
    model,
    specs,
    by.vars = NULL,
    pairwise = FALSE,
    adjust = "bonferroni"
) {
  emm.grid <- emmeans::emmeans(
    object = model,
    specs = stats::as.formula(paste0("~ ", specs)),
    by = by.vars
  )

  emm.response <- emmeans::regrid(emm.grid, transform = "response")

  if (!pairwise) {
    return(
      as_tibble(summary(emm.response, infer = c(TRUE, TRUE)))
    )
  }

  contrast.tbl <- as_tibble(
    summary(
      pairs(emm.response, adjust = adjust),
      infer = c(TRUE, TRUE)
    )
  )

  if (all(c("asymp.LCL", "asymp.UCL") %in% names(contrast.tbl))) {
    contrast.tbl <- contrast.tbl %>%
      rename(
        lower.CL = asymp.LCL,
        upper.CL = asymp.UCL
      )
  }

  return(contrast.tbl)
}
