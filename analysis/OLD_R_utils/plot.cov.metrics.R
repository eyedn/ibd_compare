###############################################################################
#           Aydin Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           plot.cov.metrics.R
###############################################################################

library(tidyverse)

plot.cov.metrics <- function(
    cov.metrics, stat, true.map, unrel.status, pops = NULL   
    ) {
  stat.sym     <- sym(stat)
  
  df <- cov.metrics %>%
    filter(
      true.map == !!true.map,
      unrel.status == !!unrel.status
    )
  
  if (!is.null(pops)) {
    df <- df %>% filter(pop %in% pops)
  }
  
  # drop rows with NA for this stat
  df <- df %>%
    filter(!is.na(!!stat.sym))
  
  if (nrow(df) == 0) {
    message("No data for stat=", stat,
            ", true.map=", true.map,
            ", unrel.status=", unrel.status,
            ". Returning NULL.")
    return(NULL)
  }
  
  # nice fixed order for infer.map if possible
  infer.order <- c("decode", "hapmapII", "1Mb1cMflat", "pyrho")
  infer.order <- intersect(infer.order, unique(df$infer.map))
  df <- df %>%
    mutate(infer.map = factor(infer.map, levels = infer.order))
  
  # y-axis label
  ylab <- switch(
    stat,
    sensitivity = "Sensitivity (TP / (TP + FN))",
    specificity = "Specificity (TN / (TN + FP))",
    PPV         = "PPV (TP / (TP + FP))",
    NPV         = "NPV (TN / (TN + FN))",
    FPR         = "FPR (FP / (FP + TN))",
    accuracy    = "Accuracy ((TP + TN) / Total)",
    stat
  )
  
  title_txt <- paste0(
    "Metric: ", stat,
    " | true.map = ", true.map,
    " | unrel.status = ", unrel.status
  )
  
  p <- ggplot(
    df,
    aes(
      x     = interaction(infer.map, pop),
      y     = !!stat.sym,
      color = pop
    )
  ) +
    # all 10 reps of same combo share same x (vertical “stack”)
    geom_point(
      position = position_jitter(width = 0.1, height = 0),
      alpha = 0.8, aes(shape = infer.map)
    ) +
    facet_wrap(~ method, nrow = 1) +
    labs(
      title = title_txt,
      y     = ylab,
      color = "population"
    ) +
    scale_color_manual(
      values = c(
        "YRI" = "#0072B2",       # blue
        "CEU" = "#D55E00",       # vermillion
        "CHB" = "#009E73"        # greenish-blue
        )
      )+
    theme_bw() +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1),
      strip.text   = element_text(size = 9)
    )
  
  p
}
