# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# run.one.confusion.analysis.task.R
# ______________________________________________________________________________


run.one.confusion.analysis.task <- function(
    metric,
    question,
    question.label,
    fixed.effects,
    contrast.var,
    by.vars,
    data,
    adj.method
) {
  message(sprintf("Starting %s %s", metric, question))

  result.tbl <- run.pairwise.question(
    data = data,
    fixed.effects = fixed.effects,
    contrast.var = contrast.var,
    by.vars = by.vars,
    metric = metric,
    question = question,
    question.label = question.label,
    adjust = adj.method
  )

  message(sprintf("Finished %s %s", metric, question))

  return(result.tbl)
}
