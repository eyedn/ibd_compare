# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# run.overall.confusion.analysis.R
# ______________________________________________________________________________


# run the full overall-only confusion analysis workflow for one metric ----
run.overall.confusion.analysis <- function(
    long.data,
    metric,
    adjust = "bonferroni",
    lrt.cores = 1
) {
  overall.data <- long.data %>%
    filter(threshold == "overall")

  overall.test.specs <- build.confusion.test.specs("overall")

  overall.fit <- run.confusion.lrt(
    data = overall.data,
    metric = metric,
    model.scope = "overall",
    test.specs = overall.test.specs,
    adjust = adjust,
    cores = lrt.cores
  )

  overall.model <- overall.fit$full.model

  emm.method <- tidy.confusion.emmeans(
    overall.model,
    "software.with.params"
  )
  emm.map <- tidy.confusion.emmeans(overall.model, "imap")
  emm.filter <- tidy.confusion.emmeans(overall.model, "filt")
  emm.pop <- tidy.confusion.emmeans(overall.model, "pop")

  method.pairs <- tidy.confusion.emmeans(
    overall.model,
    "software.with.params",
    pairwise = TRUE,
    adjust = adjust
  )
  map.pairs <- tidy.confusion.emmeans(
    overall.model,
    "imap",
    pairwise = TRUE,
    adjust = adjust
  )
  filter.pairs <- tidy.confusion.emmeans(
    overall.model,
    "filt",
    pairwise = TRUE,
    adjust = adjust
  )
  pop.pairs <- tidy.confusion.emmeans(
    overall.model,
    "pop",
    pairwise = TRUE,
    adjust = adjust
  )

  list(
    metric = metric,
    data = overall.data,
    model = overall.model,
    lrt.tbl = overall.fit$lrt.tbl,
    or.tbl = tidy.confusion.odds.ratios(overall.model),
    emm.method = emm.method,
    emm.map = emm.map,
    emm.filter = emm.filter,
    emm.pop = emm.pop,
    method.pairs = method.pairs,
    map.pairs = map.pairs,
    filter.pairs = filter.pairs,
    pop.pairs = pop.pairs
  )
}
