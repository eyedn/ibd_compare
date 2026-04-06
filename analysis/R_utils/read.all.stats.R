# ______________________________________________________________________________
# Aydin Karatas 
# ___
# University of Southern California 
# Department of Quantitative and Computational Biology 
# Mooney Lab 
# ___
# read.all.stats.R 
# ______________________________________________________________________________

library(tidyverse)

read.all.stats <- function(file, param.map) {
  all.stats <- readr::read_tsv(file, col_names = TRUE, show_col_types = FALSE
  )
  all.stats <- all.stats %>%
    mutate(
      software.code = param.map[software.with.params]
    )
  return(all.stats)
}