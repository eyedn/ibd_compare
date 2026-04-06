###############################################################################
#           plot.ibd.sharing.summary.R
###############################################################################

library(tidyverse)

plot.ibd.sharing.summary <- function(
    results.df,   # expects per-replicate results (e.g., all.res)
    true.map,
    unrel.status = c("all", "unrel"),
    value        = c("cM", "bp"),
    infer.order  = c("decode", "hapmapII", "defaultflat",
                     "1Mb1cMflat", "pyrhoYRI", "pyrhoCEU", "pyrhoCHB", "pyrho"),
    title        = NULL
) {
  unrel.status <- match.arg(unrel.status)
  value        <- match.arg(value)
  
  df <- results.df %>%
    filter(
      true.map     == !!true.map,
      unrel.status == !!unrel.status,
      value        == !!value
    ) %>%
    filter(!is.na(frob.norm))
  
  if (nrow(df) == 0) {
    message("No data for true.map = ", true.map,
            ", unrel.status = ", unrel.status,
            ", value = ", value, ". Returning NULL.")
    return(NULL)
  }
  
  infer.levels <- intersect(infer.order, unique(df$infer.map))
  df <- df %>%
    mutate(infer.map = factor(infer.map, levels = infer.levels))
  
  # --- New title formatting consistent with confusion-matrix plots ----
  # if (is.null(title)) {
  #   title <- "IBD sharing distance (Frobenius norm)"
  # }
  
  subtitle_txt <- paste0(
    "Metric: Dist. by norm | true.map = ", true.map,
    " | unrel.status = ", unrel.status,
    " | value = ", value
  )
  
  ggplot(
    df,
    aes(x = interaction(infer.map, pop), y = frob.norm, color = pop)
  ) +
    geom_point(
      position = position_jitter(width = 0.15, height = 0),
      size     = 2,
      alpha    = 0.8, aes(shape = infer.map)
    ) +
    facet_wrap(~ method, nrow = 1) +
    
    # --- Colorblind-friendly palette (Okabe–Ito) ---
    scale_color_manual(
      values = c(
        "YRI" = "#0072B2",   # blue
        "CEU" = "#D55E00",   # vermillion
        "CHB" = "#009E73"    # bluish green
      )
    ) +
    
    labs(
      title    = title,
      subtitle = subtitle_txt,
      y        = "Distance from true sharing)",
      color    = "Population"
    ) +
    theme_bw(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text  = element_text(face = "bold"),
      plot.title  = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5)
    )
}
