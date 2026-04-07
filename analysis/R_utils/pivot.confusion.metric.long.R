# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# pivot.confusion.metric.long.R
# ______________________________________________________________________________


library(tidyverse)

# reshape confusion-matrix counts into a long table for one metric ----
pivot.confusion.metric.long <- function(all.stats.per.rep, metric, thr.levels) {
  stopifnot(metric %in% c("sensitivity", "precision"))

  threshold.tbl <- tibble(
    threshold = thr.levels,
    suffix = c("", paste0(".", thr.levels[-1]))
  )

  if (metric == "sensitivity") {
    successes.prefix <- "TP"
    failures.prefix <- "FN"
    estimate.prefix <- "sensitivity"
  } else {
    successes.prefix <- "TP"
    failures.prefix <- "FP"
    estimate.prefix <- "PPV"
  }

  long.tbl <- purrr::map_dfr(
    seq_len(nrow(threshold.tbl)),
    function(idx) {
      threshold.value <- threshold.tbl$threshold[idx]
      suffix.value <- threshold.tbl$suffix[idx]

      successes.col <- paste0(successes.prefix, suffix.value)
      failures.col <- paste0(failures.prefix, suffix.value)
      estimate.col <- paste0(estimate.prefix, suffix.value)

      all.stats.per.rep %>%
        transmute(
          software = software,
          software.with.params = software.with.params,
          software.code = software.code,
          imap = imap,
          filt = filt,
          pop = pop,
          rep = rep,
          replicate.id = interaction(pop, rep, drop = TRUE),
          threshold = threshold.value,
          metric = metric,
          successes = .data[[successes.col]],
          failures = .data[[failures.col]],
          estimate = .data[[estimate.col]]
        )
    }
  )

  long.tbl %>%
    mutate(
      threshold = factor(threshold, levels = thr.levels),
      software = factor(software),
      software.with.params = factor(software.with.params),
      software.code = factor(software.code),
      imap = factor(imap),
      filt = factor(filt),
      pop = factor(pop),
      replicate.id = factor(replicate.id)
    )
}
