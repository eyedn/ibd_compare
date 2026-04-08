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
fit.confusion.glmm <- function(data, model.formula) {
  stopifnot(inherits(model.formula, "formula"))

  lme4::glmer(
    formula = model.formula,
    data = data,
    family = binomial(),
    control = lme4::glmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 2e5)
    )
  )
}
