# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology 
# Mooney Lab
# ___
# ibd_score.R
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

# read in statistics aggregated across replicates ----
all.stats <- readr::read_tsv(
  glue("{ibd.comp.dir}/true_vs_inference/all_stats.tsv"), 
  col_names = TRUE, show_col_types = FALSE
)
all.stats <- all.stats %>%
  mutate(
    software.code = param.map[software.with.params]
  )

# analyze ibd score rel true performance ----
df.reltrue <- all.stats.per.rep %>%
  filter(
    imap == "hapmap",
    filt %in% c("ibdends")
  ) %>%
  pivot_longer(
    cols = matches("^ibd\\.score\\.rel\\.true(\\.|$)"),
    values_to = "ibd.score.rel.true"
  )

ggplot(df.reltrue, aes(x = filt, y = ibd.score.rel.true, color = software.code)) +
  geom_boxplot(
    aes(group = interaction(filt, software.code)), fill = "lightgrey"
  ) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(
    x = "Filtering procedure",
    y = "Relative IBD score (inferred / true)",
    color = "Software",
    title = "Inferred vs True IBD score using ibd-ends"
  ) +
  coord_trans(y = "log10") +
  facet_wrap(~ pop) +
  theme_bw()

## plot true ibd scores ----
df.true.ibd <- all.stats.per.rep %>%
  filter(
    imap == "hapmap",
    rep == 1,
    filt == "ibdends",
    software.code == "h1"
  ) %>%
  pivot_longer(
    cols = matches("^ibd\\.score\\.err(\\.|$)"),
    values_to = "ibd.score.err"
  )
afr.score <- 119695.72
df.true.ibd <- df.true.ibd %>%
  mutate(rel.score = true.ibd.score / afr.score)
ggplot(df.true.ibd, aes(x = pop, y = true.ibd.score)) +
  geom_col(fill = "grey") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = NULL,
    y = "True IBD score",
    title = glue("IBD score"),
    fill = NULL
  ) +
  theme_bw(base_size = 32) 
ggplot(df.true.ibd, aes(x = pop, y = rel.score)) +
  geom_col(fill = "grey") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = NULL,
    y = "True IBD score rel. to AFR",
    title = glue("IBD score rel. to AFR"),
    fill = NULL
  ) +
  theme_bw(base_size = 32) 

# analyze ibd score error performance ----
df.err <- all.stats.per.rep %>%
  filter(
    imap == "hapmap",
    filt != "ibdendshot",
    # imap != "hapmap" | filt != "hot"
  ) %>%
  pivot_longer(
    cols = matches("^ibd\\.score\\.err(\\.|$)"),
    values_to = "ibd.score.err"
  )
df.err.ex <- df.err %>% filter(filt == "ibdends", imap == "hapmap", software.code %in% c("h1", "r1"))
ggplot(df.err, aes(x = pop, y = ibd.score.err, fill = software.code)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed") +
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
    x = NULL,
    y = "IBD score error",
    title = glue("IBD score error"),
    fill = NULL
  ) +
  facet_grid(filt ~ software) +
  # facet_wrap(~ filt + imap) +
  theme_bw(base_size = 32) 

# plot rel ibd score error to AFR ----
df.err <- all.stats.per.rep %>%
  # keep whatever imap/filt/software you want; don't force hapmap/ibdends if you want facets
  # filter(imap != "decode", filt != "none", imap != "hapmap" | filt != "hot") %>%
  filter(imap == "hapmap", filt == "ibdends") %>%
  group_by(imap, filt, software.code, rep) %>%   # add other ID cols here if needed (e.g., tmap)
  mutate(
    # pull the AFR baselines within each (imap,filt,software,rep) group
    true.afr = true.ibd.score[pop == "AFR"][1],
    inf.afr  = inf.ibd.score[pop == "AFR"][1],
    
    # normalize to AFR (relative-to-AFR scores)
    true.ibd.score.rel.afr = true.ibd.score / true.afr,
    inf.ibd.score.rel.afr  = inf.ibd.score  / inf.afr,
    
    # error = (true relative) - (inferred relative)
    ibd.score.rel.afr.err = true.ibd.score.rel.afr - inf.ibd.score.rel.afr
  ) %>%
  ungroup() %>%
  # optional safety: drop groups where AFR is missing or 0 (would create Inf/NA)
  filter(is.finite(true.ibd.score.rel.afr),
         is.finite(inf.ibd.score.rel.afr),
         is.finite(ibd.score.rel.afr.err))

ggplot(df.err, aes(x = pop, y = ibd.score.rel.afr.err, fill = software.code)) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(
    values = c(
      "h1" = "lightblue1",
      "h2" = "dodgerblue",
      "h3" = "blue",
      "r1" = "mistyrose",
      "r2" = "red",
      "r3" = "darkred"
    )
  ) +
  labs(
    x = NULL,
    y = "Relative IBD score error",
    title = "IBD score error (relative to AFR)" ,
    fill = NULL
  ) +
  # facet_grid(filt ~ imap) +
  # facet_wrap(~ filt + imap) +
  theme_bw(base_size = 32)
