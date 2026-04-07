# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# extract.confusion.contrasts.R
# ______________________________________________________________________________


library(tidyverse)
library(emmeans)

# extract tidy pairwise contrasts from an emmeans grid ----
extract.confusion.contrasts <- function(
    model,
    contrast.var,
    by.vars = NULL,
    metric,
    question,
    adjust = "holm"
) {
  emm.grid <- emmeans::emmeans(
    object = model,
    specs = stats::as.formula(paste0("~ ", contrast.var)),
    by = by.vars
  )

  emm.response <- emmeans::regrid(emm.grid, transform = "response")

  estimates.tbl <- as.data.frame(
    summary(emm.response, infer = c(TRUE, TRUE))
  ) %>%
    rename(level = all_of(contrast.var), level.estimate = prob)

  contrast.tbl <- as.data.frame(
    summary(
      pairs(emm.response, adjust = adjust),
      infer = c(TRUE, TRUE)
    )
  ) %>%
    tidyr::separate(
      contrast,
      into = c("level_1", "level_2"),
      sep = " - ",
      remove = FALSE
    )
  
  if (all(c("asymp.LCL", "asymp.UCL") %in% names(contrast.tbl))) {
    contrast.tbl <- contrast.tbl %>%
      rename(
        lower.CL = asymp.LCL,
        upper.CL = asymp.UCL
      )
  }
  
  if (!"lower.CL" %in% names(contrast.tbl)) {
    contrast.tbl$lower.CL <- NA_real_
  }
  
  if (!"upper.CL" %in% names(contrast.tbl)) {
    contrast.tbl$upper.CL <- NA_real_
  }

  join.vars <- by.vars

  level.1.tbl <- estimates.tbl %>%
    rename(level_1 = level, level_1.estimate = level.estimate)

  level.2.tbl <- estimates.tbl %>%
    rename(level_2 = level, level_2.estimate = level.estimate)

  if (length(join.vars) == 0) {
    contrast.tbl <- contrast.tbl %>%
      left_join(level.1.tbl, by = "level_1") %>%
      left_join(level.2.tbl, by = "level_2")
  } else {
    contrast.tbl <- contrast.tbl %>%
      left_join(level.1.tbl, by = c(join.vars, "level_1")) %>%
      left_join(level.2.tbl, by = c(join.vars, "level_2"))
  }

  contrast.tbl %>%
    mutate(
      metric = metric,
      question = question,
      contrast.variable = contrast.var,
      effect.direction = case_when(
        estimate > 0 ~ paste(level_1, ">", level_2),
        estimate < 0 ~ paste(level_1, "<", level_2),
        TRUE ~ paste(level_1, "=", level_2)
      )
    )
}
