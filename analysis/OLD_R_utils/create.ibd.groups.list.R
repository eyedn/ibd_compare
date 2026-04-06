###############################################################################
#           Aydin Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           create.ibd.groups.list.R
###############################################################################

library(tidyverse)

create.ibd.groups.list <- function(ibd.list, collapse.pyrho = TRUE) {
  # metadata from names
  label.info <- tibble::tibble(label = names(ibd.list)) %>%
    tidyr::separate(
      label,
      into   = c("method", "true.map", "infer.map", "pop", "rest"),
      sep    = "_",
      remove = FALSE
    ) %>%
    dplyr::mutate(
      # extract replicate from "10.repsub" part
      rep = str_match(rest, "^(\\d+)\\.repsub")[, 2] %>% as.integer(),
      
      # unrel status: ".repsub.unrels" vs ".repsub"
      unrel.status = if_else(
        str_detect(rest, "unrels$"), "unrel", "all"
        ),
      
      # collapse pyrhoCEU / pyrhoYRI / pyrhoCHB -> "pyrho" if requested
      true.map.group = if (collapse.pyrho) {
        if_else(
          str_detect(true.map, stringr::regex("^pyrho", ignore_case = TRUE)), 
          "pyrho", true.map
          )
        } else {
          true.map
        },
      
      # column value for inference map or truth label
      infer.map.or.truth = infer.map,
      
      # key for grouping + output list names
      group.id = stringr::str_c(true.map.group, rep, unrel.status, sep = "_"),
      
      # attach the corresponding df
      dat = ibd.list[label]
    )
  
  # split by grouping key
  group.list <- split(label.info, label.info$group.id)
  return(group.list)
}
