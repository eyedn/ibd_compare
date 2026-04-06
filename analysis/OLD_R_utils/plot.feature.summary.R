###############################################################################
#           Aydin Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           plot.feature.summary.R
###############################################################################

library(tidyverse)

plot.feature.summary <- function(
    df, pops = NULL, reps = NULL, methods = NULL, tmap,
    imaps = NULL, unrel.stat, feature, 
    min.val = NULL, max.val = NULL, bin.size = NULL, bin.by.freq = TRUE
) {
  # filter data down to specification in arguments
  if (!is.null(pops)) df <- df %>% filter(pop %in% pops)
  if (!is.null(reps)) df <- df %>% filter(rep %in% reps)
  if (!is.null(methods)) df <- df %>% filter(method %in% methods)
  if (!is.null(imaps)) df <- df %>% filter(infer.map %in% imaps)
  
  df <- df %>%filter(true.map == tmap, unrel.status == unrel.stat)
  
  # create histogram binning
  if (is.null(min.val)) min.val <- min(df[[feature]], na.rm = TRUE)
  if (is.null(max.val)) max.val <- max(df[[feature]], na.rm = TRUE)
  
  if (!is.null(bin.size)) {
    # use user-specified bin size
    if (min.val == max.val) {
      # expand a bit to avoid degenerate breaks
      min.val <- min.val - bin.size / 2
      max.val <- max.val + bin.size / 2
    }
    bin.start  <- floor(min.val / bin.size) * bin.size
    bin.end    <- ceiling(max.val / bin.size) * bin.size
    bin.breaks <- seq(bin.start, bin.end, by = bin.size)
    bin.width  <- bin.size
  } else {
    # default behavior: pretty breaks
    bin.breaks <- pretty(c(min.val, max.val), n = 30)
    bin.width  <- diff(bin.breaks)[1]
  }
  
  # group by the facet-defining variables:
  # pop, rep, infer.map, method
  group_vars <- c("pop", "rep", "infer.map", "method")
  
  hist.data <- df %>%
    group_by(across(all_of(group_vars))) %>%
    mutate(
      bin = cut(
        .data[[feature]],
        breaks = bin.breaks,
        include.lowest = TRUE,
        right = FALSE
      )
    ) %>%
    count(across(all_of(group_vars)), bin) %>%
    group_by(across(all_of(group_vars))) %>%
    mutate(freq = n / sum(n)) %>%
    ungroup() %>%
    filter(!is.na(bin)) %>%
    mutate(
      bin_low = as.numeric(sub("\\[([^,]+),.*", "\\1", as.character(bin))),
      bin_mid = bin_low + bin.width / 2
    )
  
  # facet by interaction of pop, rep, infer.map, method
  ifelse(bin.by.freq == TRUE, y.axis <- "freq", y.axis <- "n")
  histplot <- ggplot(hist.data, aes(x = bin_mid, y = .data[[y.axis]])) +
    geom_col(width = bin.width) +
    facet_wrap(
      ~ interaction(pop, rep, infer.map, method, drop = TRUE),
      scales = "fixed"
    ) +
    labs(
      x = feature,
      y = y.axis,
      title = paste0(
        "feature: ", feature,
        " | true.map = ", tmap,
        " | unrel.status = ", unrel.stat
      )
    ) +
    theme_bw() +
    theme(legend.position = "none")
  
  return(histplot)
}
