# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# fit.confusion.glmm.R
# ______________________________________________________________________________


library(lme4)

# fit a binomial glmm for confusion-matrix counts ----
fit.confusion.glmm <- function(data, fixed.effects) {
  stopifnot(is.character(fixed.effects), length(fixed.effects) == 1)

  model.formula <- as.formula(
    paste0(
      "cbind(successes, failures) ~ ",
      fixed.effects,
      " + (1 | replicate.id)"
    )
  )

  if (nrow(data) > 3000) {
    glm.formula <- as.formula(
      paste0("cbind(successes, failures) ~ ", fixed.effects)
    )

    return(
      stats::glm(
        formula = glm.formula,
        data = data,
        family = binomial()
      )
    )
  }

  fit.try <- try(
    lme4::glmer(
      formula = model.formula,
      data = data,
      family = binomial(),
      control = lme4::glmerControl(
        optimizer = "bobyqa",
        optCtrl = list(maxfun = 2e5)
      )
    ),
    silent = TRUE
  )

  if (!inherits(fit.try, "try-error")) {
    return(fit.try)
  }

  glm.formula <- as.formula(
    paste0("cbind(successes, failures) ~ ", fixed.effects)
  )

  warning(
    paste(
      "glmer failed for fixed effects",
      fixed.effects,
      "- falling back to glm for this subset."
    )
  )

  stats::glm(
    formula = glm.formula,
    data = data,
    family = binomial()
  )
}
