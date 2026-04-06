###############################################################################
#           Aydin Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           compare.ibd.sharing.all.reps.R
###############################################################################

library(tidyverse)
source("/home1/karatas/ibd_compare/analysis/true_v_inferred/make.ibd.sharing.matrix.R")
source("/home1/karatas/ibd_compare/analysis/true_v_inferred/frobenius.norm.R")
source("/home1/karatas/ibd_compare/analysis/true_v_inferred/compare.ibd.sharing.one.rep.R")

compare.ibd.sharing.all.reps <- function(
    ibd.df,
    pop,
    true.map,
    reps = 1:10,
    unrel.status = c("all", "unrel"),
    truth.method,
    truth.infer.map = NULL,
    methods,
    infer.maps,
    value = c("cM", "bp")
) {
  unrel.status <- match.arg(unrel.status)
  value        <- match.arg(value)
  
  map_dfr(
    reps,
    function(rp) {
      compare.ibd.sharing.one.rep(
        ibd.df          = ibd.df,
        pop             = pop,
        true.map        = true.map,
        rep             = rp,
        unrel.status    = unrel.status,
        truth.method    = truth.method,
        truth.infer.map = truth.infer.map,
        methods         = methods,
        infer.maps      = infer.maps,
        value           = value
      )
    }
  )
}
