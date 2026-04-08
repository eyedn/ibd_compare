# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# run.confusion.lrt.R
# ______________________________________________________________________________


library(tidyverse)

# fit full and reduced confusion models and summarize lrt results ----
run.confusion.lrt <- function(
    data,
    metric,
    model.scope,
    test.specs,
    adjust = "bonferroni"
) {
  full.terms <- test.specs$full.terms[[1]]
  full.formula <- build.confusion.model.formula(full.terms)
  full.model <- fit.confusion.glmm(data, full.formula)

  lrt.tbl <- purrr::map_dfr(
    seq_len(nrow(test.specs)),
    function(i) {
      tested.term <- test.specs$tested.term[i]
      reduced.terms <- test.specs$reduced.terms[[i]]
      reduced.formula <- build.confusion.model.formula(reduced.terms)
      reduced.model <- fit.confusion.glmm(data, reduced.formula)
      anova.tbl <- as.data.frame(
        stats::anova(reduced.model, full.model, test = "Chisq")
      )

      tibble(
        metric = metric,
        model.scope = model.scope,
        tested.term = tested.term,
        full.formula = deparse(full.formula),
        reduced.formula = deparse(reduced.formula),
        df.diff = anova.tbl$Df[2] - anova.tbl$Df[1],
        aic.diff = anova.tbl$AIC[1] - anova.tbl$AIC[2],
        lrt.statistic = anova.tbl$Chisq[2],
        p.value = anova.tbl$`Pr(>Chisq)`[2]
      )
    }
  ) %>%
    mutate(
      adjusted.p.value = p.adjust(p.value, method = adjust)
    )

  list(
    full.model = full.model,
    lrt.tbl = lrt.tbl
  )
}
