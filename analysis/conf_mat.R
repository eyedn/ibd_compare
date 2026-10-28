# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# conf_mat.R
# ______________________________________________________________________________


library(tidyverse)
library(glue)

# Resolve sources from either supported working directory.
ANALYSIS_DIR <- if (file.exists("analysis/R_utils/constants.R")) {
  "analysis"
} else "."
source(file.path(ANALYSIS_DIR, "R_utils/constants.R"), local = TRUE)
source(file.path(ANALYSIS_DIR, "R_utils/figure.core.R"), local = TRUE)
source(file.path(ANALYSIS_DIR, "R_utils/read.all.stats.R"), local = TRUE)
source(file.path(ANALYSIS_DIR, "R_utils/read.all.stats.per.rep.R"), local = TRUE)
source(file.path(ANALYSIS_DIR, "R_utils/pivot.confusion.metric.long.R"), local = TRUE)
source(file.path(ANALYSIS_DIR, "R_utils/plot.confusion.boxplots.R"), local = TRUE)

# define constants ----
per.rep.path <- glue("{ibd.comp.dir}/true_vs_inference/all_stats.per_rep.tsv")

# read in replicate-level statistics ----
all.stats.per.rep <- read.all.stats.per.rep(
  per.rep.path,
  stats.meta.names, param.map
)

validate.figure.stats(all.stats.per.rep, ANALYSIS_PROFILES$code,
                      ANALYSIS_POPULATIONS, filt.levels, FOCAL_REPLICATE)

# build long-format metric tables ----
sens.long <- pivot.confusion.metric.long(
  all.stats.per.rep = all.stats.per.rep,
  metric = "sensitivity",
  thr.levels = thr.levels
  )
precision.long <- pivot.confusion.metric.long(
  all.stats.per.rep = all.stats.per.rep,
  metric = "precision",
  thr.levels = thr.levels
  )

sens.long.overall <- sens.long %>%
  filter(threshold == "overall")
precision.long.overall <- precision.long %>%
  filter(threshold == "overall")

sens.long.bins <- sens.long %>%
  filter(threshold != "overall")
precision.long.bins <- precision.long %>%
  filter(threshold != "overall")

# summary statistics for boxplot-ready confusion metrics ----
confusion.boxplot.summary <- bind_rows(sens.long, precision.long) %>%
  group_by(pop, software.code, imap, filt, threshold, metric) %>%
  summarise(
    min = if (all(is.na(estimate))) NA_real_ else
      min(estimate, na.rm = TRUE),
    q1 = quantile(estimate, 0.25, na.rm = TRUE, names = FALSE),
    median = median(estimate, na.rm = TRUE),
    q3 = quantile(estimate, 0.75, na.rm = TRUE, names = FALSE),
    max = if (all(is.na(estimate))) NA_real_ else
      max(estimate, na.rm = TRUE),
    mean = if (all(is.na(estimate))) NA_real_ else
      mean(estimate, na.rm = TRUE),
    sd = sd(estimate, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(min, q1, median, q3, max, mean, sd),
    names_to = "statistic",
    values_to = "value"
  )


sens.all <- sens.long
precision.all <- precision.long
confusion.plots <- list()
for (pair.name in rownames(PROFILE_PAIRINGS)) {
  hapibd.code <- PROFILE_PAIRINGS[pair.name, "hapibd"]
  rapid.code <- PROFILE_PAIRINGS[pair.name, "rapid"]
  selected.codes <- c(hapibd.code, rapid.code)
  sens.long <- sens.all %>% filter(software.code %in% selected.codes)
  precision.long <- precision.all %>% filter(software.code %in% selected.codes)
  plots <- list()
  for (pop.name in c("ADX", "AFR", "EUR")) {
    for (map.name in c("hapmap", "decode", "1Mb1cM")) {
      plots[[paste("sensitivity", pop.name, map.name, sep = "_")]] <-
        plot.confusion.boxplots(sens.long, pop.name, map.name, "Sensitivity")
      plots[[paste("precision", pop.name, map.name, sep = "_")]] <-
        plot.confusion.boxplots(precision.long, pop.name, map.name, "Precision")
    }
  }
  # plots for PEQG ----
  line.sizes = 2
  filt.labels.tmp <- c(
    "none" = "None",
    "hot" = "Z",
    "ibdends" = "I",
    "hotibdends" = "Z then I",
    "ibdendshot" = "I then Z"
  )
  map.labels.tmp <- c(
    "1Mb1cM" = "flat",
    "hapmap" = "HapMapII"
  )

  # sensitivity
  sens.long.tmp <- sens.long %>%
    filter(
      software.code %in% selected.codes, threshold == "overall", pop == "ADX",
      imap != "decode", filt %in% c("none", "hot", "ibdends")
      ) %>%
    mutate(imap = factor(imap, levels = c("hapmap", "1Mb1cM")))

  sens.box <- ggplot(sens.long.tmp, aes(x = software.code, y = estimate)) +
    geom_violin(
      aes(
        group = software.code, fill = software.code,
        color = software.code
        ), linewidth = 0.75
      ) +
    geom_boxplot(
      aes(
        group = software.code, color = software.code
      ), linewidth = 0.75, width = 0.2, outliers = FALSE, fill = "white"
    ) +
    facet_grid(
      imap ~ filt,
      labeller = labeller(imap = map.labels.tmp, filt = filt.labels.tmp)
      ) +
    scale_color_manual(
      name = NULL,
      values = setNames(c("#C75DAA", "#009B9F"), selected.codes),
      labels = setNames(c("hap-IBD", "RaPID"), selected.codes)
    ) +
    scale_fill_manual(
      name = NULL,
      values = setNames(c("#EBC9DE", "#B0DDDE"), selected.codes),
      labels = setNames(c("hap-IBD", "RaPID"), selected.codes)
    ) +
    labs(
      x = "Filtering Strategy", y = "Sensitivity",
      title = NULL, fill = NULL
      ) +
    theme_bw(
      base_size = 24,
      base_line_size = line.sizes,
      base_rect_size = line.sizes
    ) +
    theme(
      panel.grid.minor = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      strip.background = element_rect(
        fill = "#fbc2c1",
        color = "black"
      ),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.box = "horizontal"
    ) +
    guides(
      fill = guide_legend(
        nrow = 1,
        byrow = TRUE
      )
    )



  # precision
  precision.long.tmp <- precision.long %>%
    filter(
      software.code %in% selected.codes, threshold == "overall", pop == "EUR",
      imap != "decode", filt %in% c("none", "hot", "ibdends")
    ) %>%
    mutate(imap = factor(imap, levels = c("hapmap", "1Mb1cM")))

  precision.box <- ggplot(
    precision.long.tmp, aes(x = software.code, y = estimate)
    ) +
    geom_violin(
      aes(
        group = software.code, fill = software.code,
        color = software.code
      ), linewidth = 0.75
    ) +
    geom_boxplot(
      aes(
        group = software.code, color = software.code
      ), linewidth = 0.75, width = 0.2, outliers = FALSE, fill = "white"
    ) +
    facet_grid(
      imap ~ filt,
      labeller = labeller(imap = map.labels.tmp, filt = filt.labels.tmp)
    ) +
    scale_color_manual(
      name = NULL,
      values = setNames(c("#C75DAA", "#009B9F"), selected.codes),
      labels = setNames(c("hap-IBD", "RaPID"), selected.codes)
    ) +
    scale_fill_manual(
      name = NULL,
      values = setNames(c("#EBC9DE", "#B0DDDE"), selected.codes),
      labels = setNames(c("hap-IBD", "RaPID"), selected.codes)
    ) +
    scale_y_continuous(limits = c(0,1), breaks = c(0, .25, 0.5, 0.75, 1)) +
    labs(
      x = "Filtering Strategy", y = "Precision",
      title = NULL, fill = NULL
    ) +
    theme_bw(
      base_size = 24,
      base_line_size = line.sizes,
      base_rect_size = line.sizes
    ) +
    theme(
      panel.grid.minor = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      strip.background = element_rect(
        fill = "#fbc2c1",
        color = "black"
      ),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.box = "horizontal"
    ) +
    guides(
      fill = guide_legend(
        nrow = 1,
        byrow = TRUE
      )
    )



  # sensitivity and precision but with different levels
  thr.labels.tmp <- c(
    "overall" = "overall",
    "2"   = "2–2.5",
    "2.5" = "2.5–3",
    "3"   = "3–3.5"
  )
  metric.labels.tmp <- c(
    "sensitivity" = "TPR",
    "precision" = "PPV"
  )
  pd <- position_dodge(width = 0.8)

  thresh.long <- rbind(sens.long, precision.long) %>%
    filter(
      software.code %in% selected.codes,
      threshold %in% c("overall", "2", "2.5", "3"),
      pop == "ADX",
      imap == "hapmap",
      filt %in% c("hot", "ibdends")
    ) %>%
    mutate(
      imap = factor(imap, levels = c("hapmap", "1Mb1cM")),
      metric = factor(metric, levels = c("sensitivity", "precision"))
        )

  thresh.box <- ggplot(
    thresh.long,
    aes(
      x = factor(threshold),
      y = estimate,
      fill = software.code,
      color = software.code
    )
  ) +
    geom_violin(
      position = pd,
      linewidth = 0.75
    ) +
    geom_boxplot(
      position = pd,
      width = 0.2,
      linewidth = 0.75,
      fill = "white",
      outlier.shape = NA
    ) +
    facet_grid(
      metric ~ filt,
      labeller = labeller(metric = metric.labels.tmp, filt = filt.labels.tmp),
      scales = "free_y"
    ) +
    scale_x_discrete(labels = thr.labels.tmp) +
    scale_color_manual(
      name = NULL,
      values = setNames(c("#C75DAA", "#009B9F"), selected.codes),
      labels = setNames(c("hap-IBD", "RaPID"), selected.codes)
    ) +
    scale_fill_manual(
      name = NULL,
      values = setNames(c("#EBC9DE", "#B0DDDE"), selected.codes),
      labels = setNames(c("hap-IBD", "RaPID"), selected.codes)
    ) +
    labs(
      x = "IBD Length bins (cM)", y = "Value with HapMapII",
      title = NULL, fill = NULL
    ) +
    theme_bw(
      base_size = 24,
      base_line_size = line.sizes,
      base_rect_size = line.sizes
    ) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 28, hjust = 1, vjust = 1),
      axis.title.x = element_blank(),
      # axis.text.x = element_blank(),
      # axis.ticks.x = element_blank(),
      strip.background = element_rect(
        fill = "#fbc2c1",
        color = "black"
      ),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.box = "horizontal"
    ) +
    guides(
      fill = guide_legend(
        nrow = 1,
        byrow = TRUE
      )
    )



  plots$sensitivity.presentation <- sens.box
  plots$precision.presentation <- precision.box
  plots$threshold.presentation <- thresh.box
  confusion.plots[[pair.name]] <- plots
}
confusion.plots <- organize.figure.plots(confusion.plots)
invisible(confusion.plots)
