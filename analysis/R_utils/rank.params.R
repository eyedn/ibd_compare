# ______________________________________________________________________________
# Aydin Karatas 
# ___
# University of Southern California 
# Department of Quantitative and Computational Biology 
# Mooney Lab 
# ___
# rank.params.R 
# ______________________________________________________________________________

library(dplyr)

rank.params <- function(df, metric, direction = c("max", "min")) {
  direction <- match.arg(direction)
  
  df %>%
    group_by(pop) %>%
    arrange(
      if (direction == "max") desc(.data[[metric]]) else .data[[metric]],
      .by_group = TRUE
    ) %>%
    mutate(rank = row_number()) %>%
    ungroup()
}
