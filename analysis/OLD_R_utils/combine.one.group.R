###############################################################################
#           Aydin Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           combine.one.group.R
###############################################################################

library(tidyverse)

combine.one.group <- function(group.df) {
  purrr::map2_dfr(
    group.df$dat,
    seq_len(nrow(group.df)),
    ~ dplyr::mutate(
      .x,
      infer.map.or.truth = group.df$infer.map.or.truth[.y],
      pop                = group.df$pop[.y]
    )
  )
}