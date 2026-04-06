###############################################################################
#           Aydin Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           compare.ibd.sharing.one.rep.R
###############################################################################

library(tidyverse)
source("/home1/karatas/ibd_compare/analysis/true_v_inferred/make.ibd.sharing.matrix.R")
source("/home1/karatas/ibd_compare/analysis/true_v_inferred/frobenius.norm.R")

compare.ibd.sharing.one.rep <- function(
    ibd.df,
    pop,
    true.map,
    rep,
    unrel.status = c("all", "unrel"),
    truth.method,
    truth.infer.map = NULL,
    methods,
    infer.maps,
    value = c("cM", "bp")
) {
  unrel.status <- match.arg(unrel.status)
  value        <- match.arg(value)
  
  # For true calls, infer.map is the same as the truth.method (e.g. "true2Mb")
  if (is.null(truth.infer.map)) {
    truth.infer.map <- truth.method
  }
  
  ## -------------------- build truth sharing matrix ---------------------------
  truth.res <- make.ibd.sharing.matrix(
    ibd.df       = ibd.df,
    pop          = pop,
    true.map     = true.map,
    rep          = rep,
    unrel.status = unrel.status,
    method       = truth.method,
    infer.map    = truth.infer.map,
    value        = value
  )
  
  truth.mat   <- truth.res$matrix
  truth.samps <- truth.res$samples
  
  # if no truth samples / matrix, bail out cleanly
  if (length(truth.samps) == 0 || nrow(truth.mat) == 0) {
    return(tibble())
  }
  
  ## -------------------- compare each (method, infer.map) combo --------------
  combos <- expand_grid(method = methods, infer.map = infer.maps)
  
  res <- pmap_dfr(
    combos,
    function(method, infer.map) {
      
      # for pyrho maps, return early if the map does not match the pop
      if (str_detect(infer.map, "pyrho") && !str_detect(infer.map, pop)) {
        return(tibble())
      }
      
      infer.res <- make.ibd.sharing.matrix(
        ibd.df       = ibd.df,
        pop          = pop,
        true.map     = true.map,
        rep          = rep,
        unrel.status = unrel.status,
        method       = method,
        infer.map    = infer.map,
        value        = value,
        samples      = truth.samps
      )
      
      infer.mat <- infer.res$matrix
      
      # collapse pyrho<POP> -> "pyrho" in the *output*, like compute.confusion.metrics
      infer.map.out <- if (str_detect(infer.map, "pyrho")) "pyrho" else infer.map
      
      if (!is.matrix(infer.mat) ||
          any(dim(infer.mat) != dim(truth.mat))) {
        return(
          tibble(
            pop          = pop,
            true.map     = true.map,
            rep          = rep,
            unrel.status = unrel.status,
            method       = method,
            infer.map.out,
            value        = value,
            frob.norm    = NA_real_
          ) %>%
            rename(infer.map = infer.map.out)
        )
      }
      
      tibble(
        pop          = pop,
        true.map     = true.map,
        rep          = rep,
        unrel.status = unrel.status,
        method       = method,
        infer.map.out,
        value        = value,
        frob.norm    = frobenius.norm(truth.mat, infer.mat)
      ) %>%
        rename(infer.map = infer.map.out)
    }
  )
  
  res
}
