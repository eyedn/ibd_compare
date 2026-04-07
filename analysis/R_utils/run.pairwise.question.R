# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# run.pairwise.question.R
# ______________________________________________________________________________


library(tidyverse)

# fit one question-level model and extract pairwise contrasts ----
run.pairwise.question <- function(
    data,
    fixed.effects,
    contrast.var,
    by.vars = NULL,
    metric,
    question,
    question.label,
    adjust = "holm"
) {
  model <- fit.confusion.glmm(data, fixed.effects)

  extract.confusion.contrasts(
    model = model,
    contrast.var = contrast.var,
    by.vars = by.vars,
    metric = metric,
    question = question,
    adjust = adjust
  ) %>%
    mutate(question.label = question.label) %>%
    complete.condition.cols()
}
