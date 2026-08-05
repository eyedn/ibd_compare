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
library(lme4)
library(emmeans)

# source functions ----
source("R_utils/constants.R")
source("R_utils/log.msg.R")
source("R_utils/read.all.stats.per.rep.R")
source("R_utils/pivot.confusion.metric.long.R")
source("R_utils/build.confusion.model.formula.R")
source("R_utils/build.confusion.test.specs.R")
source("R_utils/fit.confusion.glmm.R")
source("R_utils/run.confusion.lrt.R")
source("R_utils/tidy.confusion.emmeans.R")
source("R_utils/tidy.confusion.odds.ratios.R")
source("R_utils/run.overall.confusion.analysis.R")
source("R_utils/run.bin.effects.confusion.analysis.R")
source("R_utils/plot.confusion.boxplots.R")

# define constants ----
adj.method <- "bonferroni"
per.rep.path <- glue("{ibd.comp.dir}/true_vs_inference/all_stats.per_rep.tsv")

# read in replicate-level statistics ----
all.stats.per.rep <- read.all.stats.per.rep(
  per.rep.path,
  stats.meta.names, param.map
)

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
    min = min(estimate, na.rm = TRUE),
    q1 = quantile(estimate, 0.25, na.rm = TRUE, names = FALSE),
    median = median(estimate, na.rm = TRUE),
    q3 = quantile(estimate, 0.75, na.rm = TRUE, names = FALSE),
    max = max(estimate, na.rm = TRUE),
    mean = mean(estimate, na.rm = TRUE),
    sd = sd(estimate, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(min, q1, median, q3, max, mean, sd),
    names_to = "statistic",
    values_to = "value"
  )

# fit four primary glmm models and omnibus lrt tables ----
# note: format of data from pivot.confusion.metric.long
# sens.analysis <- run.overall.confusion.analysis(
#   long.data = sens.long.overall, 
#   adjust.method = adj.method
#   )
# 
# precision.analysis <- run.overall.confusion.analysis(
#   long.data = precision.long.overall, 
#   adjust.method = adj.method
#   )
# 
# sens.bin.analysis <- run.bin.effects.confusion.analysis(
#   long.data = sens.long.bins,
#   adjust.method = adj.method
#   )
# 
# precision.bin.analysis <- run.bin.effects.confusion.analysis(
#   long.data = precision.long.bins,
#   adjust.method = adj.method
#   )

# build manuscript-facing plotting tables ----
sens.long.tmp <- sens.long %>% filter(software.code %in% c("h3", "r1"), threshold %in% c("overall", "2"))
plot.confusion.boxplots(sens.long, "ADX", "hapmap", "Sensitivity")
plot.confusion.boxplots(sens.long, "AFR", "hapmap", "Sensitivity")
# plot.confusion.boxplots(sens.long, "ASIA", "hapmap", "Sensitivity")
plot.confusion.boxplots(sens.long, "EUR", "hapmap", "Sensitivity")

plot.confusion.boxplots(sens.long, "ADX", "decode", "Sensitivity")
plot.confusion.boxplots(sens.long, "AFR", "decode", "Sensitivity")
# plot.confusion.boxplots(sens.long, "ASIA", "decode", "Sensitivity")
plot.confusion.boxplots(sens.long, "EUR", "decode", "Sensitivity")

plot.confusion.boxplots(sens.long, "ADX", "1Mb1cM", "Sensitivity")
plot.confusion.boxplots(sens.long, "AFR", "1Mb1cM", "Sensitivity")
# plot.confusion.boxplots(sens.long, "ASIA", "1Mb1cM", "Sensitivity")
plot.confusion.boxplots(sens.long, "EUR", "1Mb1cM", "Sensitivity")

plot.confusion.boxplots(precision.long, "ADX", "hapmap", "Precision")
plot.confusion.boxplots(precision.long, "AFR", "hapmap", "Precision")
# plot.confusion.boxplots(precision.long, "ASIA", "hapmap", "Precision")
plot.confusion.boxplots(precision.long, "EUR", "hapmap", "Precision")

plot.confusion.boxplots(precision.long, "ADX", "decode", "Precision")
plot.confusion.boxplots(precision.long, "AFR", "decode", "Precision")
# plot.confusion.boxplots(precision.long, "ASIA", "decode", "Precision")
plot.confusion.boxplots(precision.long, "EUR", "decode", "Precision")

plot.confusion.boxplots(precision.long, "ADX", "1Mb1cM", "Precision")
plot.confusion.boxplots(precision.long, "AFR", "1Mb1cM", "Precision")
# plot.confusion.boxplots(precision.long, "ASIA", "1Mb1cM", "Precision")
plot.confusion.boxplots(precision.long, "EUR", "1Mb1cM", "Precision")

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
    software.code %in% c("h3", "r1"), threshold == "overall", pop == "ADX",
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
    values = c("h3" = "#C75DAA", "r1" = "#009B9F"),
    labels = c("h3" = "hap-IBD", "r1" = "RaPID")
  ) +
  scale_fill_manual(
    name = NULL,
    values = c("h3" = "#EBC9DE", "r1" = "#B0DDDE"),
    labels = c("h3" = "hap-IBD", "r1" = "RaPID")
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

if (interactive()) {
  print(sens.box)
}

pdf(
  "~/proj/ibd_compare_results/sens_box.pdf",
  width = 6,
  height = 5,
  family = "Helvetica"
)

print(sens.box)

dev.off()

# precision
precision.long.tmp <- precision.long %>% 
  filter(
    software.code %in% c("h3", "r1"), threshold == "overall", pop == "EUR", 
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
    values = c("h3" = "#C75DAA", "r1" = "#009B9F"),
    labels = c("h3" = "hap-IBD", "r1" = "RaPID")
  ) +
  scale_fill_manual(
    name = NULL,
    values = c("h3" = "#EBC9DE", "r1" = "#B0DDDE"),
    labels = c("h3" = "hap-IBD", "r1" = "RaPID")
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

if (interactive()) {
  print(precision.box)
}

pdf(
  "~/proj/ibd_compare_results/precision_box.pdf",
  width = 6,
  height = 5,
  family = "Helvetica"
)

print(precision.box)

dev.off()

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
    software.code %in% c("h3", "r1"), 
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
    values = c("h3" = "#C75DAA", "r1" = "#009B9F"),
    labels = c("h3" = "hap-IBD", "r1" = "RaPID")
  ) +
  scale_fill_manual(
    name = NULL,
    values = c("h3" = "#EBC9DE", "r1" = "#B0DDDE"),
    labels = c("h3" = "hap-IBD", "r1" = "RaPID")
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

if (interactive()) {
  print(thresh.box)
}

pdf(
  "~/proj/ibd_compare_results/thresh_box.pdf",
  width = 6,
  height = 5,
  family = "Helvetica"
)

print(thresh.box)

dev.off()
