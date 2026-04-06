# ______________________________________________________________________________
# Aydin Karatas
# ---
# University of Southern California
# Department of Quantitative and Computational Biology 
# Mooney Lab
# ___
# subset.replicate.R
###############################################################################

library(tidyverse)

subset.replicate <- function(ibd.fname, rep.fname, ibd.calls, pop.ids) {
  # determine which pop these ibd calls are from
  pop <- names(pop.ids)[
    sapply(names(pop.ids), function(name) grepl(name, ibd.fname))
  ]
  
  # read in replicate subset and isolate samples form this pop
  subset <- readr::read_csv(
    file = rep.fname, col_names = TRUE, show_col_types = FALSE
  )
  subset <- subset %>% 
    rename(indiv.id = individual_id, pop.id = population_id) %>%
    filter(pop.id == pop.ids[pop])
  
  # apply pop-specific, replicate subset to ibd calls 
  ibd.calls <- ibd.calls %>% filter(
    sam1 %in% subset$indiv.id & sam2 %in% subset$indiv.id 
    & sam1 != sam2
  )
  print(paste0("replicate subset filtered calls: ", nrow(ibd.calls)))
  
  return(ibd.calls)
}
