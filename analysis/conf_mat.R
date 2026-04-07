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
library(furrr)
library(future)
library(progressr)
library(glue)
library(stats)
library(emmeans)

# source function ----
source("R_utils/constants.R")
source("R_utils/read.all.stats.R")
source("R_utils/read.all.stats.per.rep.R")

# read in statistics aggregated across reps / per rep ----
all.stats <- read.all.stats(
  glue("{ibd.comp.dir}/true_vs_inference/all_stats.tsv"), param.map
  )
all.stats.per.rep <- read.all.stats.per.rep(
  glue("{ibd.comp.dir}/true_vs_inference/all_stats.per_rep.tsv"), 
  stats.meta.names, param.map
  )

# define focal variables ----
focal.pop <- "ADMIX"
focal.imap <- "hapmap"
focal.software <- "hapibd"

# analyze sensitivity ----
## pop specific sensitivity ----
df.sens <- all.stats.per.rep %>%
  filter(
    pop == focal.pop,
    imap == focal.imap,
    filt != "ibdendshot"
  ) %>%
  pivot_longer(
    cols = matches("^sensitivity(\\.|$)"),
    names_to = "threshold",
    values_to = "sensitivity"
  ) %>%
  mutate(
    threshold = case_when(
      threshold == "sensitivity" ~ "overall",
      TRUE ~ sub("^sensitivity\\.", "", threshold)
    ),
    threshold = factor(threshold, levels = thr.levels)
  )

ggplot(df.sens, aes(x = threshold, y = sensitivity)) +
  geom_boxplot(
    aes(group = interaction(threshold, software.code), fill = software.code)
  ) +
  facet_grid(filt ~ software) +
  scale_x_discrete(labels = thr.labels) +
  scale_fill_manual(
    values = c(
      "h1" = "lightblue1",
      "h2"  = "dodgerblue",
      "h3" = "blue",
      "r1" = "mistyrose",
      "r2"  = "red",
      "r3" = "darkred"
    )
  ) +
  labs(
    x = "IBD length bins (cM)",
    y = "Sensitivity",
    title = glue("Sensitivty of {focal.pop} IBD with a {focal.imap} map"),
    fill = NULL
  ) +
  theme_bw(base_size = 32) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

## across pops sensitivity ----
df.sens.allpops <- all.stats.per.rep %>%
  filter(
    # software == focal.software,
    software.code == "h1",
    imap == focal.imap,
    # filt != "ibdendshot"
    filt == "ibdends"
  ) %>%
  pivot_longer(
    cols = matches("^sensitivity(\\.|$)"),
    names_to = "threshold",
    values_to = "sensitivity"
  ) %>%
  mutate(
    threshold = case_when(
      threshold == "sensitivity" ~ "overall",
      TRUE ~ sub("^sensitivity\\.", "", threshold)
    ),
    threshold = factor(threshold, levels = thr.levels)
  )

ggplot(df.sens.allpops, aes(x = threshold, y = sensitivity)) +
  geom_boxplot(
    aes(group = interaction(threshold, software.code), fill = software.code)
  ) +
  facet_grid(filt ~ pop) +
  scale_x_discrete(labels = thr.labels) +
  scale_fill_manual(
    values = c(
      "h1" = "lightblue1",
      "h2"  = "dodgerblue",
      "h3" = "blue",
      "r1" = "mistyrose",
      "r2"  = "red",
      "r3" = "darkred"
    )
  ) +
  labs(
    x = "IBD length bins (cM)",
    y = "Sensitivity",
    title = glue("Sensitivty of {focal.software} calls with a {focal.imap} map"),
    fill = NULL
  ) +
  theme_bw(base_size = 32) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

## for gsg symposium ----
df.sens.allpops.subthresh <- df.sens.allpops %>% 
  filter(threshold %in% c("overall", "2"))
ggplot(
  df.sens.allpops.subthresh, aes(x = threshold, y = sensitivity)
  ) +
  geom_boxplot(
    aes(group = interaction(threshold, pop), fill = pop), 
    outliers = FALSE,
    linewidth = 1
  ) +
  scale_x_discrete(labels = thr.labels.descriptive) +
  scale_fill_manual(
    values = c(
      "ADMIX"  = "dodgerblue",
      "AFR" = "mistyrose",
      "ASIA"  = "red",
      "EUR" = "darkred"
    )
  ) +
  labs(
    y = "True positive rate",
    title = "Simulated IBD Inference performance",
    fill = "Populations"
  ) +
  theme_bw(base_size = 42, base_line_size = 2, base_rect_size = 3) +
  theme(
    text = element_text(family = "Arial", color = "black"),
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title.x = element_blank(),
    legend.position = c(0.05, 0.05),
    legend.justification = c(0, 0),
    legend.background = element_rect(
      fill = alpha("white"), color = "black", linewidth = 2
      )
  )

# analyze ppv ----
## pop specific ppv ----
df.ppv <- all.stats.per.rep %>%
  filter(
    pop == focal.pop,
    imap == focal.imap,
    filt != "ibdendshot"
  ) %>%
  pivot_longer(
    cols = matches("^ppv(\\.|$)"),
    names_to = "threshold",
    values_to = "ppv"
  ) %>%
  mutate(
    threshold = case_when(
      threshold == "ppv" ~ "overall",
      TRUE ~ sub("^ppv\\.", "", threshold)
    ),
    threshold = factor(threshold, levels = thr.levels)
  )

ggplot(df.ppv, aes(x = threshold, y = ppv)) +
  geom_boxplot(
    aes(group = interaction(threshold, software.code), fill = software.code)
  ) +
  facet_grid(filt ~ software) +
  scale_x_discrete(labels = thr.labels) +
  scale_fill_manual(
    values = c(
      "h1" = "lightblue1",
      "h2"  = "dodgerblue",
      "h3" = "blue",
      "r1" = "mistyrose",
      "r2"  = "red",
      "r3" = "darkred"
    )
  ) +
  labs(
    x = "IBD length bins (cM)",
    y = "ppv",
    title = glue("ppv of {focal.pop} IBD with a {focal.imap} map"),
    # title = glue("ppv of hap-IBD + ibd-ends across populations"),
    fill = NULL
  ) +
  theme_bw(base_size = 32) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

## across pops ppv ----
df.ppv.allpops <- all.stats.per.rep %>%
  filter(
    software == focal.software,
    imap == focal.imap,
    filt != "ibdendshot"
  ) %>%
  pivot_longer(
    cols = matches("^ppv(\\.|$)"),
    names_to = "threshold",
    values_to = "ppv"
  ) %>%
  mutate(
    threshold = case_when(
      threshold == "ppv" ~ "overall",
      TRUE ~ sub("^ppv\\.", "", threshold)
    ),
    threshold = factor(threshold, levels = thr.levels)
  )

ggplot(df.ppv.allpops, aes(x = threshold, y = ppv)) +
  geom_boxplot(
    aes(group = interaction(threshold, software.code), fill = software.code)
  ) +
  facet_grid(filt ~ pop) +
  scale_x_discrete(labels = thr.labels) +
  scale_fill_manual(
    values = c(
      "h1" = "lightblue1",
      "h2"  = "dodgerblue",
      "h3" = "blue",
      "r1" = "mistyrose",
      "r2"  = "red",
      "r3" = "darkred"
    )
  ) +
  labs(
    x = "IBD length bins (cM)",
    y = "ppv",
    title = glue("ppv of {focal.software} calls with a {focal.imap} map"),
    fill = NULL
  ) +
  theme_bw(base_size = 32) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
