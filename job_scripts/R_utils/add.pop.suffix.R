# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology 
# Mooney Lab
# ___
# add.pop.suffix.R
# ______________________________________________________________________________

library(tidyverse)
library(tools)

add.pop.suffix <- function(path, pop) {
  
  # If Pyrho path → REMOVE "_pop"
  if (grepl("pyrho", path, ignore.case = TRUE)) {
    # remove the literal `_pop` substring
    new.path <- sub("_pop", "", path, fixed = TRUE)
    return(new.path)
  }
  
  # Otherwise → replace `_pop` with `_<pop>`
  new.path <- sub("_pop", paste0("_", pop), path, fixed = TRUE)
  
  return(new.path)
}
