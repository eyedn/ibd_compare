# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology 
# Mooney Lab
# ___
# read.all.stats.per.rep.R
# ______________________________________________________________________________


library(tidyverse)

read.all.stats.per.rep <- function(file, meta.names, param.map) {
  all.stats.per.rep <- readr::read_tsv(
    file, col_names = TRUE, show_col_types = FALSE
  )
  stopifnot(ncol(all.stats.per.rep) >= length(meta.names))
  all.stats.per.rep <- all.stats.per.rep %>%
    mutate(
      software.code = param.map[software.with.params]
    )
  return(all.stats.per.rep)
}