###############################################################################
# compare.ibd.sharing.all.pops.maps.R
###############################################################################
library(tidyverse)
library(furrr)

source("/home1/karatas/ibd_compare/analysis/true_v_inferred/compare.ibd.sharing.all.reps.R")

compare.ibd.sharing.all.pops.maps <- function(
    ibd.df,
    pops,
    true.maps,
    reps         = 1:10,
    unrel.status = c("all", "unrel"),
    truth.method,
    truth.infer.map = NULL,  # if NULL, compare.ibd.sharing.one.rep() will default to truth.method
    methods,
    infer.maps,
    value        = c("cM", "bp"),
    .progress    = TRUE
) {
  unrel.status <- match.arg(unrel.status)
  value        <- match.arg(value)
  
  # grid of (pop, true.map)
  combo.grid <- expand_grid(
    pop      = pops,
    true.map = true.maps
  )
  
  # parallel map over (pop, true.map)
  all.res <- future_pmap_dfr(
    combo.grid,
    function(pop, true.map) {
      compare.ibd.sharing.all.reps(
        ibd.df          = ibd.df,
        pop             = pop,
        true.map        = true.map,
        reps            = reps,
        unrel.status    = unrel.status,
        truth.method    = truth.method,
        truth.infer.map = truth.infer.map,
        methods         = methods,
        infer.maps      = infer.maps,
        value           = value
      )
    },
    .progress = .progress
  )
  
  all.res
}
