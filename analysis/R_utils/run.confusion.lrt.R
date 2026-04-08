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
library(future)
library(furrr)
library(progressr)

# fit full and reduced confusion models and summarize lrt results ----
run.confusion.lrt <- function(
    data,
    metric,
    model.scope,
    test.specs,
    adjust = "bonferroni",
    cores = 1
) {
  full.terms <- test.specs$full.terms[[1]]
  full.formula <- build.confusion.model.formula(full.terms)
  full.model <- fit.confusion.glmm(data, full.formula)

  reduced.specs <- tibble(
    tested.term = test.specs$tested.term,
    reduced.terms = test.specs$reduced.terms
  ) %>%
    mutate(
      reduced.formula = purrr::map(
        reduced.terms,
        build.confusion.model.formula
      )
    )

  future::plan(future::multisession, workers = cores)
  on.exit(future::plan(future::sequential), add = TRUE)

  reduced.fit.tbl <- progressr::with_progress({
    p <- progressr::progressor(along = seq_len(nrow(reduced.specs)))

    reduced.models <- furrr::future_map(
      seq_len(nrow(reduced.specs)),
      function(i) {
        tested.term <- reduced.specs$tested.term[i]
        reduced.formula <- reduced.specs$reduced.formula[[i]]

        message(sprintf("Starting reduced model for %s", tested.term))
        reduced.model <- fit.confusion.glmm(data, reduced.formula)
        message(sprintf("Finished reduced model for %s", tested.term))
        p(sprintf("Completed reduced model for %s", tested.term))

        reduced.model
      },
      .options = furrr::furrr_options(seed = TRUE)
    )

    reduced.specs %>%
      mutate(
        reduced.model = reduced.models
      )
  })

  lrt.tbl <- purrr::map_dfr(
    seq_len(nrow(reduced.fit.tbl)),
    function(i) {
      tested.term <- reduced.fit.tbl$tested.term[i]
      reduced.formula <- reduced.fit.tbl$reduced.formula[[i]]
      reduced.model <- reduced.fit.tbl$reduced.model[[i]]
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
