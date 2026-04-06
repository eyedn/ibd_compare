# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology 
# Mooney Lab
# ___
# combine.ibd.list.R
# ______________________________________________________________________________

library(tidyverse)

combine.ibd.list <- function(ibd.list, collapse.pyrho = TRUE) {
  
  meta <- tibble(label = names(ibd.list)) %>%
    # split label into base label + params
    separate(
      label,
      into   = c("base_label", "params"),
      sep    = "-",
      remove = FALSE,
      fill   = "right",
      extra  = "merge"   # keep all remaining '-' in params if any
    ) %>%
    mutate(
      params = replace_na(params, "")
    ) %>%
    # parse the base label like before
    separate(
      base_label,
      into   = c("method", "true.map", "infer.map", "pop", "rest"),
      sep    = "_",
      remove = FALSE,
      fill   = "right",
      extra  = "merge"
    ) %>%
    mutate(
      # replicate index if present 
      rep = suppressWarnings(
        as.integer(
          if_else(
            str_detect(rest %||% "", "\\."),
            str_extract(rest %||% "", "^[^.]+"),
            rest %||% ""
          )
        )
      ),
      
      # unrel status
      unrel.status = if_else(
        str_detect(rest %||% "", "unrels$"), "unrel", "all"
      ),
      
      # collapse pyrhoPOP into pyrho
      true.map.group = if (collapse.pyrho) {
        if_else(
          str_detect(true.map %||% "", regex("^pyrho", ignore_case = TRUE)),
          "pyrho",
          true.map
        )
      } else {
        true.map
      }
    )
  
  combined <- purrr::imap_dfr(
    ibd.list,
    ~ {
      info <- dplyr::filter(meta, label == .y) %>% dplyr::slice(1)
      
      .x %>%
        mutate(
          method       = info$method,
          pop          = info$pop,
          rep          = info$rep,
          unrel.status = info$unrel.status,
          true.map     = info$true.map.group,
          infer.map    = info$infer.map,
          params        = info$params
        )
    }
  )
  
  combined
}
