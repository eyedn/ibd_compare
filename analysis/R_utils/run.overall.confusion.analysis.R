# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# run.overall.confusion.analysis.R
# ______________________________________________________________________________

library(lme4)
library(tidyr)

# run the full overall-only confusion analysis workflow for sens/precision ----
run.overall.confusion.analysis <- function(
    long.data,
    adjust.method = "bonferroni"
) {
  # note: format of data from pivot.confusion.metric.long
  # fit binomial glmm
  log.msg("fitting model")
  model <- lme4::glmer(
    cbind(successes, failures) ~
      software.with.params + imap + filt + pop + (1 | replicate.id),
    data = long.data,
    family = binomial(),
    control = lme4::glmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 2e5)
    )
  )
  
  # run omnibus likelihood ratio tests for each factor
  log.msg("running omnibus lrt for each factor")
  mod.no.method <- update(model, . ~ . - software.with.params)
  mod.no.map <- update(model, . ~ . - imap)
  mod.no.filter <- update(model, . ~ . - filt)
  mod.no.pop <- update(model, . ~ . - pop)
  lrt.method <- anova(mod.no.method, model, test = "Chisq")
  lrt.map <- anova(mod.no.map, model, test = "Chisq")
  lrt.filter <- anova(mod.no.filter, model, test = "Chisq")
  lrt.pop <- anova(mod.no.pop, model, test = "Chisq")
  
  # collect omnibus lrt into one table
  log.msg("generating lrt table")
  lrt.tbl <- tibble::tibble(
    term = c("software.with.params", "imap", "filt", "pop"),
    lrt.statistic = c(
      lrt.method$Chisq[2],
      lrt.map$Chisq[2],
      lrt.filter$Chisq[2],
      lrt.pop$Chisq[2]
    ),
    p.value = c(
      lrt.method$`Pr(>Chisq)`[2],
      lrt.map$`Pr(>Chisq)`[2],
      lrt.filter$`Pr(>Chisq)`[2],
      lrt.pop$`Pr(>Chisq)`[2]
    )
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
  
  # adjusted post hoc comparisons for any significant factor
  log.msg("generating methods emmeans")
  emm.method <- emmeans::emmeans(
    model,
    ~ software.with.params,
    type = "response"
  )
  emm.method.pairs <- pairs(emm.method, adjust = adjust.method)
  summary.emm.method.pairs <- summary(emm.method.pairs)
  summary.emm.method.pairs <- summary.emm.method.pairs %>% 
    mutate(log.odds.ratio = log(odds.ratio))
    
  log.msg("generating maps emmeans")
  emm.map <- emmeans::emmeans(
    model,
    ~ imap,
    type = "response"
  )
  emm.map.pairs <- pairs(emm.map, adjust = adjust.method)
  summary.emm.map.pairs <- summary(emm.map.pairs)
  summary.emm.map.pairs <- summary.emm.map.pairs %>% 
    mutate(log.odds.ratio = log(odds.ratio))
  
  log.msg("generating filter emmeans")
  emm.filter <- emmeans::emmeans(
    model,
    ~ filt,
    type = "response"
  )
  emm.filter.pairs <- pairs(emm.filter, adjust = adjust.method)
  summary.emm.filter.pairs <- summary(emm.filter.pairs)
  summary.emm.filter.pairs <- summary.emm.filter.pairs %>% 
    mutate(log.odds.ratio = log(odds.ratio))
  
  log.msg("generating pop emmeans")
  emm.pop <- emmeans::emmeans(
    model,
    ~ pop,
    type = "response"
  )
  emm.pop.pairs <- pairs(emm.pop, adjust = adjust.method)
  summary.emm.pop.pairs <- summary(emm.pop.pairs)
  summary.emm.pop.pairs <- summary.emm.pop.pairs %>% 
    mutate(log.odds.ratio = log(odds.ratio))
  
  # return analysis results
  return(
    list(
      "lrt.tbl" = lrt.tbl,
      "or.tbl" = or.tbl,
      "summary.emm.method.pairs" = summary.emm.method.pairs,
      "summary.emm.map.pairs" = summary.emm.map.pairs,
      "summary.emm.filter.pairs" = summary.emm.filter.pairs,
      "summary.emm.pop.pairs" = summary.emm.pop.pairs
      )
    )
}
