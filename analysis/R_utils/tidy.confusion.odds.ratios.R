# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# tidy.confusion.odds.ratios.R
# ______________________________________________________________________________


library(tidyverse)
library(lme4)

# extract fixed-effect odds ratios and wald confidence intervals ----
tidy.confusion.odds.ratios <- function(model) {
  tibble(
    term = names(lme4::fixef(model)),
    estimate = lme4::fixef(model)
  ) %>%
    mutate(
      odds.ratio = exp(estimate)
    ) %>%
    bind_cols(
      as.data.frame(confint(model, parm = "beta_", method = "Wald"))
    ) %>%
    rename(
      ci.lo.logit = `2.5 %`,
      ci.hi.logit = `97.5 %`
    ) %>%
    mutate(
      ci.lo.or = exp(ci.lo.logit),
      ci.hi.or = exp(ci.hi.logit)
    )
}
