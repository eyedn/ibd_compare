###############################################################################
#           Aydin Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           combine.ibd.cov.list.R
###############################################################################

library(tidyverse)

combine.ibd.cov.list <- function(ibd.list) {
  meta <- tibble(label = names(ibd.list)) %>%
    separate(
      label,
      into   = c(
        "X1", "X2", "true.map", "X3", "pop", "rep", "X4", "method", "X5", 
        "infer.map", "X6", "X7"
      ),
      sep    = "_",
      remove = FALSE
    ) %>%
    select(-X1, -X2, -X3, -X4, -X5, -X6, -X7) %>%
    mutate(
      # unrel status
      unrel.status = if_else(
        str_detect(rep, "unrels$"), "unrel", "all"
      ),
      
      # replicate index from "10.repsub*"
      rep = str_match(rep, "^(\\d+)\\.repsub")[, 2] %>% as.integer(),
    )
  
  # combine all dfs
  combined <- purrr::imap_dfr(
    ibd.list,
    ~ {
      info <- filter(meta, label == .y) %>% slice(1)
      
      .x %>%
        mutate(
          method     = info$method,          
          pop          = info$pop,
          rep          = info$rep,
          unrel.status = info$unrel.status,
          true.map     = info$true.map,
          infer.map    = info$infer.map
        )
    }
  )
  
  return(combined)
}
