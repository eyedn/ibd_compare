# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology 
# Mooney Lab
# ___
# subset.unrels.R
# ______________________________________________________________________________


library(tidyverse)

subset.unrels <- function(unrel.fname, ibd.calls,
                          sample.id.mode = c("simulation", "empirical")) {
  sample.id.mode <- match.arg(sample.id.mode)
  message("filtering ibd calls for unrelated in ", unrel.fname)
  
  # read in unrelateds
  unrels <- readr::read_tsv(
    file = unrel.fname, 
    col_names = FALSE, show_col_types = FALSE
  )
  unrels <- unrels %>%
    rename(indiv.id = X1) %>%
    mutate(indiv.id = as.character(indiv.id))

  if (sample.id.mode == "simulation") {
    unrels <- unrels %>%
      mutate(indiv.id = str_extract(indiv.id, "\\d+$"))
  }

  ibd.calls <- ibd.calls %>%
    mutate(
      sam1 = as.character(sam1),
      sam2 = as.character(sam2)
    )
  
  # subset for unrels
  ibd.calls <- ibd.calls %>% filter(
    sam1 %in% unrels$indiv.id & sam2 %in% unrels$indiv.id 
    & sam1 != sam2
  )
  message("unrels subset filtered calls: ", nrow(ibd.calls))
  
  return(list("calls" = ibd.calls, "n.unrels" = nrow(unrels)))
}
