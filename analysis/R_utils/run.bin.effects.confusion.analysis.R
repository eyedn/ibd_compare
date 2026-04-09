# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# run.bin.effects.confusion.analysis.R
# ______________________________________________________________________________


library(lme4)
library(tidyr)

# run the bin-effects confusion analysis workflow for sens/precision ----
run.bin.effects.confusion.analysis <- function(
    long.data,
    adjust.method = "bonferroni"
) {
  # fit binomial glmm with bin as the main fixed effect of interest
  log.msg("fitting model")
  model <- lme4::glmer(
    cbind(successes, failures) ~ bin + method + map + filter + population +
      (1 | replicate),
    data = long.data,
    family = binomial(),
    control = lme4::glmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 2e5)
    )
  )

  # test whether bin is significant after adjusting for covariates
  log.msg("running omnibus lrt for bin")
  mod.no.bin <- update(model, . ~ . - bin)
  lrt.bin <- anova(mod.no.bin, model, test = "Chisq")

  # collect the omnibus lrt for bin
  log.msg("generating lrt table")
  lrt.tbl <- tibble::tibble(
    term = "bin",
    lrt.statistic = lrt.bin$Chisq[2],
    p.value = lrt.bin$`Pr(>Chisq)`[2]
  ) %>%
    dplyr::mutate(
      adjusted.p.value = p.adjust(p.value, method = adjust.method)
    )

  # get odds ratios and Wald CIs for the fixed effects
  log.msg("generating or table")
  or.tbl <- tibble::tibble(
    term = names(lme4::fixef(model)),
    estimate = lme4::fixef(model)
  ) %>%
    dplyr::mutate(
      odds.ratio = exp(estimate)
    ) %>%
    dplyr::bind_cols(
      as.data.frame(
        confint(model, parm = "beta_", method = "Wald")
      )
    ) %>%
    dplyr::rename(
      log.odds.rato = estimate,
      ci.lo.logit = `2.5 %`,
      ci.hi.logit = `97.5 %`
    ) %>%
    dplyr::mutate(
      ci.lo.or = exp(ci.lo.logit),
      ci.hi.or = exp(ci.hi.logit)
    ) %>%
    dplyr::select(term, log.odds.rato, odds.ratio, ci.lo.or, ci.hi.or)

  # adjusted post hoc comparisons for bin levels
  log.msg("generating bin emmeans")
  emm.bin <- emmeans::emmeans(
    model,
    ~ bin,
    type = "response"
  )
  emm.bin.pairs <- pairs(emm.bin, adjust = adjust.method)
  summary.emm.bin.pairs <- summary(emm.bin.pairs) %>%
    dplyr::mutate(log.odds.ratio = log(odds.ratio))

  # return analysis results
  return(
    list(
      "lrt.tbl" = lrt.tbl,
      "or.tbl" = or.tbl,
      "summary.emm.bin.pairs" = summary.emm.bin.pairs
    )
  )
}
