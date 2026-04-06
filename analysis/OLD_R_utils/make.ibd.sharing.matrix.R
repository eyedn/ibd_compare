###############################################################################
#           Aydin Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           make.ibd.sharing.matrix.R
###############################################################################

library(tidyverse)

make.ibd.sharing.matrix <- function(
    ibd.df,
    pop,
    true.map,
    rep,
    unrel.status = c("all", "unrel"),
    method,
    infer.map,
    value        = c("cM", "bp"),
    samples      = NULL,
    denom.bp     = 200e6   # approximate genome length per pair, in bp
) {
  unrel.status <- match.arg(unrel.status)
  value        <- match.arg(value)
  
  len.col <- if (value == "bp") "bp.len" else "cM.len"
  
  df.sub <- ibd.df %>%
    filter(
      pop          == !!pop,
      true.map     == !!true.map,
      rep          == !!rep,
      unrel.status == !!unrel.status,
      method       == !!method,
      infer.map    == !!infer.map
    )
  
  if (is.null(samples)) {
    samples <- sort(unique(c(df.sub$sam1, df.sub$sam2)))
  }
  samples <- as.integer(samples)
  n <- length(samples)
  
  if (n == 0) {
    return(list(matrix = matrix(numeric(0), 0, 0), samples = integer(0)))
  }
  
  if (nrow(df.sub) == 0) {
    M <- matrix(0, nrow = n, ncol = n,
                dimnames = list(samples, samples))
    return(list(matrix = M, samples = samples))
  }
  
  pairs <- df.sub %>%
    mutate(
      id1 = pmin(sam1, sam2),
      id2 = pmax(sam1, sam2)
    ) %>%
    group_by(id1, id2) %>%
    summarise(
      total = sum(.data[[len.col]], na.rm = TRUE),
      .groups = "drop"
    )
  
  M <- matrix(0, nrow = n, ncol = n,
              dimnames = list(samples, samples))
  
  i <- match(pairs$id1, samples)
  j <- match(pairs$id2, samples)
  
  for (k in seq_along(i)) {
    if (is.na(i[k]) || is.na(j[k])) next
    M[i[k], j[k]] <- pairs$total[k]
    M[j[k], i[k]] <- pairs$total[k]
  }
  
  # convert total length to approximate proportion of genome shared
  if (value == "bp") {
    M <- M / denom.bp
  }
  # if value == "cM", you could instead divide by a cM-based denominator
  
  list(matrix = M, samples = samples)
}
