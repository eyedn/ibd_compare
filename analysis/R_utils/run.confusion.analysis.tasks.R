# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# run.confusion.analysis.tasks.R
# ______________________________________________________________________________


library(tidyverse)
library(future)
library(furrr)
library(progressr)

run.confusion.analysis.tasks <- function(
    analysis.task.tbl,
    adj.method,
    n.workers
) {
  future::plan(future::multisession, workers = n.workers)
  on.exit(future::plan(future::sequential), add = TRUE)

  analysis.result.list <- progressr::with_progress({
    p <- progressr::progressor(along = seq_len(nrow(analysis.task.tbl)))

    furrr::future_pmap(
      .l = analysis.task.tbl,
      .f = function(
          metric,
          question,
          question.label,
          fixed.effects,
          contrast.var,
          by.vars,
          data
      ) {
        p(sprintf("%s %s", metric, question))

        run.pairwise.question(
          data = data,
          fixed.effects = fixed.effects,
          contrast.var = contrast.var,
          by.vars = by.vars,
          metric = metric,
          question = question,
          question.label = question.label,
          adjust = adj.method
        )
      },
      .options = furrr::furrr_options(seed = TRUE)
    )
  })

  bind_rows(analysis.result.list) %>%
    group_by(metric, question) %>%
    mutate(
      adjusted.p.value = p.adjust(contrast.p.value, method = adj.method)
    ) %>%
    ungroup()
}
