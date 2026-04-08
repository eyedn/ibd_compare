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
source("R_utils/read.all.stats.per.rep.R")
source("R_utils/pivot.confusion.metric.long.R")
source("R_utils/build.confusion.model.formula.R")
source("R_utils/build.confusion.test.specs.R")
source("R_utils/fit.confusion.glmm.R")
source("R_utils/run.confusion.lrt.R")
source("R_utils/tidy.confusion.emmeans.R")

# define constants ----
adj.method <- "bonferroni"
per.rep.path <- glue("{ibd.comp.dir}/true_vs_inference/all_stats.per_rep.tsv")

# read in replicate-level statistics ----
all.stats.per.rep <- read.all.stats.per.rep(
  per.rep.path,
  stats.meta.names, param.map
)

software.palette <- c(
  "h1" = "lightblue1",
  "h2" = "dodgerblue",
  "h3" = "blue",
  "r1" = "mistyrose",
  "r2" = "red",
  "r3" = "darkred"
)

pop.palette <- c(
  "ADMIX" = "dodgerblue",
  "AFR" = "mistyrose",
  "ASIA" = "red",
  "EUR" = "darkred"
)

focal.pop <- "ADMIX"
focal.imap <- "hapmap"
focal.software.code <- "h1"
focal.filter <- "ibdends"

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

# build data subsets for overall and bin models ----
overall.sens.data <- sens.long %>%
  filter(threshold == "overall")

overall.precision.data <- precision.long %>%
  filter(threshold == "overall")

bin.sens.data <- sens.long %>%
  filter(threshold != "overall")

bin.precision.data <- precision.long %>%
  filter(threshold != "overall")

# fit four primary glmm models and omnibus lrt tables ----
overall.test.specs <- build.confusion.test.specs("overall")
bin.test.specs <- build.confusion.test.specs("bin")

overall.sens.fit <- run.confusion.lrt(
  data = overall.sens.data,
  metric = "sensitivity",
  model.scope = "overall",
  test.specs = overall.test.specs,
  adjust = adj.method
)

overall.precision.fit <- run.confusion.lrt(
  data = overall.precision.data,
  metric = "precision",
  model.scope = "overall",
  test.specs = overall.test.specs,
  adjust = adj.method
)

bin.sens.fit <- run.confusion.lrt(
  data = bin.sens.data,
  metric = "sensitivity",
  model.scope = "bin",
  test.specs = bin.test.specs,
  adjust = adj.method
)

bin.precision.fit <- run.confusion.lrt(
  data = bin.precision.data,
  metric = "precision",
  model.scope = "bin",
  test.specs = bin.test.specs,
  adjust = adj.method
)

overall.sens.model <- overall.sens.fit$full.model
overall.precision.model <- overall.precision.fit$full.model
bin.sens.model <- bin.sens.fit$full.model
bin.precision.model <- bin.precision.fit$full.model

overall.sens.lrt.tbl <- overall.sens.fit$lrt.tbl
overall.precision.lrt.tbl <- overall.precision.fit$lrt.tbl
bin.sens.lrt.tbl <- bin.sens.fit$lrt.tbl
bin.precision.lrt.tbl <- bin.precision.fit$lrt.tbl

confusion.glmm.lrt.results <- bind_rows(
  overall.sens.lrt.tbl,
  overall.precision.lrt.tbl,
  bin.sens.lrt.tbl,
  bin.precision.lrt.tbl
)

# response-scale marginal means for interpretation ----
overall.sens.emmeans <- list(
  method = tidy.confusion.emmeans(overall.sens.model, "software.with.params"),
  map = tidy.confusion.emmeans(overall.sens.model, "imap"),
  filter = tidy.confusion.emmeans(overall.sens.model, "filt"),
  population = tidy.confusion.emmeans(overall.sens.model, "pop")
)

overall.precision.emmeans <- list(
  method = tidy.confusion.emmeans(
    overall.precision.model, "software.with.params"
  ),
  map = tidy.confusion.emmeans(overall.precision.model, "imap"),
  filter = tidy.confusion.emmeans(overall.precision.model, "filt"),
  population = tidy.confusion.emmeans(overall.precision.model, "pop")
)

bin.sens.emmeans <- list(
  threshold = tidy.confusion.emmeans(bin.sens.model, "threshold"),
  threshold.by.method = tidy.confusion.emmeans(
    bin.sens.model, "threshold", by.vars = c("software.with.params")
  )
)

bin.precision.emmeans <- list(
  threshold = tidy.confusion.emmeans(bin.precision.model, "threshold"),
  threshold.by.method = tidy.confusion.emmeans(
    bin.precision.model, "threshold", by.vars = c("software.with.params")
  )
)

# build manuscript-facing plotting tables ----
df.sens <- sens.long %>%
  filter(
    pop == focal.pop,
    imap == focal.imap,
    filt != "ibdendshot"
  )

ggplot(df.sens, aes(x = threshold, y = estimate)) +
  geom_boxplot(
    aes(group = interaction(threshold, software.code), fill = software.code)
  ) +
  facet_grid(filt ~ software) +
  scale_x_discrete(labels = thr.labels) +
  scale_fill_manual(values = software.palette) +
  labs(
    x = "IBD length bins (cM)",
    y = "Sensitivity",
    title = glue("Sensitivty of {focal.pop} IBD with a {focal.imap} map"),
    fill = NULL
  ) +
  theme_bw(base_size = 32) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

## across pops sensitivity ----
df.sens.allpops <- sens.long %>%
  filter(
    software.code == focal.software.code,
    imap == focal.imap,
    filt == focal.filter
  )

ggplot(df.sens.allpops, aes(x = threshold, y = estimate)) +
  geom_boxplot(
    aes(group = interaction(threshold, software.code), fill = software.code)
  ) +
  facet_grid(filt ~ pop) +
  scale_x_discrete(labels = thr.labels) +
  scale_fill_manual(values = software.palette) +
  labs(
    x = "IBD length bins (cM)",
    y = "Sensitivity",
    title = glue(
      "Sensitivty of {focal.software.code} calls with a {focal.imap} map"
    ),
    fill = NULL
  ) +
  theme_bw(base_size = 32) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

## for gsg symposium ----
# df.sens.allpops.subthresh <- df.sens.allpops %>% 
#   filter(threshold %in% c("overall", "2"))
# ggplot(
#   df.sens.allpops.subthresh, aes(x = threshold, y = estimate)
#   ) +
#   geom_boxplot(
#     aes(group = interaction(threshold, pop), fill = pop), 
#     outliers = FALSE,
#     linewidth = 1
#   ) +
#   scale_x_discrete(labels = thr.labels.descriptive) +
#   scale_fill_manual(values = pop.palette) +
#   labs(
#     y = "True positive rate",
#     title = "Simulated IBD Inference performance",
#     fill = "Populations"
#   ) +
#   theme_bw(base_size = 42, base_line_size = 2, base_rect_size = 3) +
#   theme(
#     text = element_text(family = "Arial", color = "black"),
#     axis.text.x = element_text(color = "black"),
#     axis.text.y = element_text(color = "black"),
#     axis.title.x = element_blank(),
#     legend.position = c(0.05, 0.05),
#     legend.justification = c(0, 0),
#     legend.background = element_rect(
#       fill = alpha("white"), color = "black", linewidth = 2
#       )
#   )

# analyze precision ----
## pop specific precision ----
df.ppv <- precision.long %>%
  filter(
    pop == focal.pop,
    imap == focal.imap,
    filt != "ibdendshot"
  )

ggplot(df.ppv, aes(x = threshold, y = estimate)) +
  geom_boxplot(
    aes(group = interaction(threshold, software.code), fill = software.code)
  ) +
  facet_grid(filt ~ software) +
  scale_x_discrete(labels = thr.labels) +
  scale_fill_manual(values = software.palette) +
  labs(
    x = "IBD length bins (cM)",
    y = "Precision",
    title = glue("Precision of {focal.pop} IBD with a {focal.imap} map"),
    fill = NULL
  ) +
  theme_bw(base_size = 32) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

## across pops precision ----
df.ppv.allpops <- precision.long %>%
  filter(
    software.code == focal.software.code,
    imap == focal.imap,
    filt != "ibdendshot"
  )

ggplot(df.ppv.allpops, aes(x = threshold, y = estimate)) +
  geom_boxplot(
    aes(group = interaction(threshold, software.code), fill = software.code)
  ) +
  facet_grid(filt ~ pop) +
  scale_x_discrete(labels = thr.labels) +
  scale_fill_manual(values = software.palette) +
  labs(
    x = "IBD length bins (cM)",
    y = "Precision",
    title = glue(
      "Precision of {focal.software.code} calls with a {focal.imap} map"
    ),
    fill = NULL
  ) +
  theme_bw(base_size = 32) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
