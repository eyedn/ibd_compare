# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# build.confusion.model.formula.R
# ______________________________________________________________________________


# plot boxplot of sens/prec given a pop-map specificaiton ----
plot.confusion.boxplots <- function(
    data.long, focal.pop, focal.imap, metric.name
    ) {
  df2plot <- data.long %>%
    filter(
      pop == focal.pop,
      imap == focal.imap
    )
  
  ggplot(df2plot, aes(x = threshold, y = estimate)) +
    geom_boxplot(
      aes(group = interaction(threshold, software.code), fill = software.code)
    ) +
    facet_grid(filt ~ .) +
    scale_x_discrete(labels = thr.labels) +
    scale_fill_manual(values = software.palette) +
    labs(
      x = "IBD length bins (cM)",
      y = metric.name,
      title = glue("{metric.name} of {focal.pop} IBD with a {focal.imap} map"),
      fill = NULL
    ) +
    theme_bw(base_size = 32) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
}
