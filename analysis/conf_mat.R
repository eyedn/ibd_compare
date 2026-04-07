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
library(future)
library(furrr)
library(glue)
library(lme4)
library(emmeans)

# source functions ----
source("R_utils/constants.R")
source("R_utils/read.all.stats.per.rep.R")
source("R_utils/pivot.confusion.metric.long.R")
source("R_utils/fit.confusion.glmm.R")
source("R_utils/extract.confusion.contrasts.R")
source("R_utils/complete.condition.cols.R")
source("R_utils/run.pairwise.question.R")

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
n.workers <- min(8, future::availableCores())

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

metric.long.list <- list(
  sensitivity = sens.long,
  precision = precision.long
)

# build question-specific tasks ----
analysis.task.tbl <- tibble(
  metric = c(
    "sensitivity", "sensitivity", "sensitivity", "sensitivity",
    "sensitivity", "precision", "precision", "precision",
    "precision", "precision"
  ),
  question = c(
    "question_1", "question_2", "question_3", "question_4", "question_5",
    "question_1", "question_2", "question_3", "question_4", "question_5"
  ),
  question.label = c(
    "Different strategies at hapmap while holding filter constant",
    "Map changes within strategy while holding filter constant",
    "Filter changes within strategy with hapmap fixed",
    "Population differences within fixed inference settings",
    "IBD bin differences within fixed inference settings",
    "Different strategies at hapmap while holding filter constant",
    "Map changes within strategy while holding filter constant",
    "Filter changes within strategy with hapmap fixed",
    "Population differences within fixed inference settings",
    "IBD bin differences within fixed inference settings"
  ),
  fixed.effects = c(
    "software.with.params * filt + pop",
    "imap * software.with.params * filt + pop",
    "filt * software.with.params + pop",
    "pop * software.with.params * imap * filt",
    "threshold * (software.with.params + imap + filt + pop)",
    "software.with.params * filt + pop",
    "imap * software.with.params * filt + pop",
    "filt * software.with.params + pop",
    "pop * software.with.params * imap * filt",
    "threshold * (software.with.params + imap + filt + pop)"
  ),
  contrast.var = c(
    "software.with.params", "imap", "filt", "pop", "threshold",
    "software.with.params", "imap", "filt", "pop", "threshold"
  ),
  by.vars = list(
    c("filt"),
    c("software.with.params", "filt"),
    c("software.with.params"),
    c("software.with.params", "imap", "filt"),
    c("software.with.params", "imap", "filt", "pop"),
    c("filt"),
    c("software.with.params", "filt"),
    c("software.with.params"),
    c("software.with.params", "imap", "filt"),
    c("software.with.params", "imap", "filt", "pop")
  ),
  data = list(
    sens.long %>% filter(threshold == "overall", imap == "hapmap"),
    sens.long %>% filter(threshold == "overall"),
    sens.long %>% filter(threshold == "overall", imap == "hapmap"),
    sens.long %>% filter(threshold == "overall"),
    sens.long %>% filter(threshold != "overall"),
    precision.long %>% filter(threshold == "overall", imap == "hapmap"),
    precision.long %>% filter(threshold == "overall"),
    precision.long %>% filter(threshold == "overall", imap == "hapmap"),
    precision.long %>% filter(threshold == "overall"),
    precision.long %>% filter(threshold != "overall")
  )
)

# run the five inferential analyses in parallel ----
plan(multisession, workers = n.workers)
analysis.result.list <- furrr::future_pmap(
  .l = analysis.task.tbl,
  .f = function(
      metric,
      question,
      question.label,
      fixed.effects,
      contrast.var,
      by.vars,
      data
  ) {
    run.pairwise.question(
      data = data,
      fixed.effects = fixed.effects,
      contrast.var = contrast.var,
      by.vars = by.vars,
      metric = metric,
      question = question,
      question.label = question.label,
      adjust = adj.method
    )
  },
  .options = furrr::furrr_options(seed = TRUE)
)
plan(sequential)

analysis.results <- bind_rows(analysis.result.list) %>%
  group_by(metric, question) %>%
  mutate(
    adjusted.p.value = p.adjust(contrast.p.value, method = adj.method)
  ) %>%
  ungroup()

confusion.metric.glmm.results <- analysis.results %>%
  arrange(metric, question, filt, software.with.params, imap, pop, threshold)

sensitivity.glmm.results <- confusion.metric.glmm.results %>%
  filter(metric == "sensitivity")

precision.glmm.results <- confusion.metric.glmm.results %>%
  filter(metric == "precision")

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
