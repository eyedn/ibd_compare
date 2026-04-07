# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology 
# Mooney Lab
# ___
# distribution.R
# ______________________________________________________________________________


library(tidyverse)
library(furrr)
library(future)
library(progressr)
library(glue)
library(stats)
library(emmeans)

# source function ----
source("R_utils/run.ks.R")
source("R_utils/constants.R")
source("R_utils/read.pop.calls.R")
source("R_utils/read.rep.calls.R")

# set up for parallel processing ----
cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
handlers("cli")

# read in unrels ----
pop.info.dir <- glue("{sim.dir}/pop_info")
unrels.tbl <- tibble(
  file = list.files(
    pop.info.dir, full.names = TRUE, 
    pattern = "^HapMapII_GRCh38_.*_unrel_samples\\.txt$"
  )
) %>%
  mutate(
    rep = as.integer(str_match(basename(file), "^HapMapII_GRCh38_(\\d+)_")[,2]),
    pop = str_match(basename(file), "^HapMapII_GRCh38_\\d+_([^_]+)_unrel_samples\\.txt$")[,2],
    ids = map(file, ~ readr::read_tsv(.x, col_names = FALSE, show_col_types = FALSE) %>%
                transmute(ids = as.integer(str_remove(X1, "^tsk_"))) %>%
                pull(ids))
  ) %>%
  select(rep, pop, ids)

# read in calls data for a specific pop ----
len.focal.pop <- "AFR"
imap.all.calls <- read.pop.calls(
  len.focal.pop, dir.paths.map, unrels.tbl, cores
  )

# build table of files for rep 1 ----
len.focal.rep <- 1
rep1.calls <- read.rep.calls(
  len.focal.rep, dir.paths.map, unrels.tbl, cores
  )

# analyze length distribution using ks testing ----
## ks testing procedure ----
ks.res <- imap.all.calls %>%   
  filter(software.code != "true", is.unrel == TRUE) %>%
  group_by(pop, rep, software.code, param.code, imap, filt) %>% 
  summarise(inf = list(cM.len), .groups = "drop") %>%
  left_join(
    imap.all.calls %>%
      filter(software.code == "true") %>%
      group_by(pop, rep) %>%
      summarise(tru = list(cM.len), .groups = "drop"),
    by = c("pop", "rep")
  ) %>%
  mutate(
    p.value = map2_dbl(tru, inf, run.ks),
    n.inf   = map_int(inf, length),
    n.true  = map_int(tru, length)
  ) %>%
  select(-tru, -inf)

ks.res <- ks.res %>%
  group_by(pop) %>%
  filter(filt != "ibdendshot") %>%
  mutate(
    p.adj = p.adjust(p.value, method = "BH"),
    log.p.adj = log10(p.adj),
    is.sig = p.adj <= 0.05,
    filt = factor(filt, levels = c("none", "hot", "ibdends")),
    imap = factor(imap, levels = c("hapmap", "decode", "1Mb1cM")),
    software.code = factor(software.code),
    param.code = factor(param.code),
  ) %>%
  ungroup()

ks.bar <- ks.res %>%
  mutate(outcome = if_else(p.adj < 0.05, "reject", "don't reject")) %>%
  count(pop, imap, filt, software.code, param.code, outcome, name = "n") %>%
  group_by(pop, imap, filt, software.code) %>%
  mutate(frac = n / sum(n)) %>%
  ungroup() %>%
  mutate(
    sw.family = case_when(
      str_starts(software.code, "h") ~ "hapibd",
      str_starts(software.code, "r") ~ "rapid",
      TRUE ~ NA_character_
    ),
    fill.key = interaction(sw.family, outcome, sep = "_")
  )

## plot fraction of null hypo. rejections ----
ks.bar.ex <- ks.bar %>% filter(filt == "ibdends", imap == "hapmap")
ggplot(
  ks.bar,
  aes(x = param.code, y = frac, fill = fill.key)
) +
  geom_col(position = "stack", width = 0.8, color = "black", linewidth = 0.5) +
  facet_grid(imap ~ filt) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_manual(
    values = c(
      "hapibd_don't reject" = "lightblue1",
      "hapibd_reject"       = "darkblue",
      "rapid_don't reject"  = "mistyrose",
      "rapid_reject"        = "darkred"
    ),
    labels = c(
      "hapibd_don't reject" = "hap-IBD: don't reject",
      "hapibd_reject"       = "hap-IBD: reject",
      "rapid_don't reject"  = "RaPID: don't reject",
      "rapid_reject"        = "RaPID: reject"
    )
  ) +
  labs(
    title = glue(
      "K-S Test True vs Inferred {len.focal.pop} IBD length distributions (unrelated pairs)"
    ),
    x = "Software / parameterization",
    y = "Fraction of replicates",
    fill = NULL
  ) +
  theme_bw(base_size = 32)

# plot ibd distributions of rep 1 ----
## plot ibd length distribution ----
sub.call.y.max <- 300 # specific to each pop
sub.calls <- rep1.calls %>%
  filter(
    param.code == "true" | 
      (param.code == "r1" & filt == "ibdends" & imap == "hapmap") | 
      (param.code == "r1" & filt == "ibdends" & imap == "1Mb1cM"),
    pop == len.focal.pop,
    is.unrel == TRUE
  ) %>%
  # hard coded facet wrap order
  mutate(
    facet = interaction(param.code, imap, filt, sep = ", "),
    facet = fct_relevel(
      facet,
      "h1, 1Mb1cM, ibdends",
      "h1, hapmap, ibdends",
      "true, true, none",
      "r1, 1Mb1cM, ibdends",
      "r1, hapmap, ibdends"
    )
  )

ggplot(sub.calls) +
  geom_histogram(
    aes(x = cM.len, fill = param.code), color = "black",  linewidth = 0.5
  ) + 
  facet_wrap(~ facet, drop = TRUE) +
  coord_cartesian(ylim = c(0, sub.call.y.max)) +
  scale_fill_manual(
    values = c(
      "h1" = "lightblue1",
      "h2"  = "dodgerblue",
      "h3" = "blue",
      "r1" = "mistyrose",
      "r2"  = "red",
      "r3" = "darkred"
    ),
    labels = NULL
  ) + 
  labs(
    title = glue("Length of {len.focal.pop} IBD segments"),
    x = "Length (cM)",
    y = "count",
    fill = NULL
  ) +
  theme_bw(base_size = 32) +
  theme(legend.position = "none")

## plot ibd colored by centromere location ----
sub.calls.cen <- sub.calls %>%
  mutate(
    not.hits.cen = (bp.end < cen.end.bp) | (bp.str > cen.end.bp),
    hits.cen = if_else(not.hits.cen, "doesn't overlap", "overlaps centromere"),
    fill.plot = if_else(hits.cen == "overlaps centromere", "cen", param.code),
    fill.plot = factor(fill.plot, levels = c("true", "h1", "r1", "cen"))
  )

ggplot(sub.calls.cen) +
  geom_histogram(
    aes(x = cM.len, fill = fill.plot),
    position = "stack",
    color = "black",
    linewidth = 0.5
  ) +
  facet_wrap(~ facet, drop = TRUE) +
  coord_cartesian(ylim = c(0, sub.call.y.max)) +
  scale_fill_manual(
    values = c(
      "h1"   = "lightblue1",
      "r1"   = "mistyrose",
      "true" = "lightgrey",
      "cen"  = "black"
    )
  ) +
  labs(
    title = glue(
      "{sub.call.focal.pop} IBD segments overlapping centromere (black bars)"
    ),
    x = "Length (cM)",
    y = "Count",
    fill = NULL
  ) +
  theme_bw(base_size = 32) +
  theme(legend.position = "none")

## plot ibd colored by before centromere ----
sub.calls.bfcen <- sub.calls %>%
  mutate(
    bef.cen = bp.end < cen.start.bp,
    bef.cen = if_else(bef.cen, "before centromere", "in/after centromere"),
    fill.plot = if_else(bef.cen == "before centromere", "bfcen", param.code),
    fill.plot = factor(fill.plot, levels = c("true", "h1", "r1", "bfcen"))
  )

ggplot(sub.calls.bfcen) +
  geom_histogram(
    aes(x = cM.len, fill = fill.plot),
    position = "stack",
    color = "black",
    linewidth = 0.5
  ) +
  facet_wrap(~ facet, drop = TRUE) +
  coord_cartesian(ylim = c(0, sub.call.y.max)) +
  scale_fill_manual(
    values = c(
      "h1"   = "lightblue1",
      "r1"   = "mistyrose",
      "true" = "lightgrey",
      "bfcen"  = "black"
    )
  ) +
  labs(
    title = glue(
      "{sub.call.focal.pop} IBD segments before centromere (black bars)"
    ),
    x = "Length (cM)",
    y = "Count",
    fill = NULL
  ) +
  theme_bw(base_size = 32) +
  theme(legend.position = "none")

## plot ibd colored by after centromere ----
sub.calls.aftcen <- sub.calls %>%
  mutate(
    aft.cen = bp.str > cen.start.bp,
    aft.cen = if_else(aft.cen, "after centromere", "before/in centromere"),
    fill.plot = if_else(aft.cen == "after centromere", "aftcen", param.code),
    fill.plot = factor(fill.plot, levels = c("true", "h1", "r1", "aftcen"))
  )

ggplot(sub.calls.aftcen) +
  geom_histogram(
    aes(x = cM.len, fill = fill.plot),
    position = "stack",
    color = "black",
    linewidth = 0.5
  ) +
  facet_wrap(~ facet, drop = TRUE) +
  coord_cartesian(ylim = c(0, sub.call.y.max)) +
  scale_fill_manual(
    values = c(
      "h1"   = "lightblue1",
      "r1"   = "mistyrose",
      "true" = "lightgrey",
      "aftcen"  = "black"
    )
  ) +
  labs(
    title = glue(
      "{len.focal.pop} IBD segments after centromere (black bars)"
    ),
    x = "Length (cM)",
    y = "Count",
    fill = NULL
  ) +
  theme_bw(base_size = 32) +
  theme(legend.position = "none")

# analyze length distribution post-centromere using ks testing ----
## ks testing procedure ----
ks.res.aftcen <- imap.all.calls %>%   
  filter(software.code != "true", is.unrel == TRUE, bp.str > cen.end.bp) %>%
  group_by(pop, rep, software.code, param.code, imap, filt) %>% 
  summarise(inf = list(cM.len), .groups = "drop") %>%
  left_join(
    imap.all.calls %>%
      filter(software.code == "true") %>%
      group_by(pop, rep) %>%
      summarise(tru = list(cM.len), .groups = "drop"),
    by = c("pop", "rep")
  ) %>%
  mutate(
    p.value = map2_dbl(tru, inf, run.ks),
    n.inf   = map_int(inf, length),
    n.true  = map_int(tru, length)
  ) %>%
  select(-tru, -inf)

ks.res.aftcen <- ks.res.aftcen %>%
  group_by(pop) %>%
  filter(filt != "ibdendshot") %>%
  mutate(
    p.adj = p.adjust(p.value, method = "BH"),
    log.p.adj = log10(p.adj),
    is.sig = p.adj <= 0.05,
    filt = factor(filt, levels = c("none", "hot", "ibdends")),
    imap = factor(imap, levels = c("hapmap", "decode", "1Mb1cM")),
    software.code = factor(software.code),
    param.code = factor(param.code),
  ) %>%
  ungroup()

ks.bar.aftcen <- ks.res.aftcen %>%
  mutate(outcome = if_else(p.adj < 0.05, "reject", "don't reject")) %>%
  count(pop, imap, filt, software.code, param.code, outcome, name = "n") %>%
  group_by(pop, imap, filt, software.code) %>%
  mutate(frac = n / sum(n)) %>%
  ungroup() %>%
  mutate(
    sw.family = case_when(
      str_starts(software.code, "h") ~ "hapibd",
      str_starts(software.code, "r") ~ "rapid",
      TRUE ~ NA_character_
    ),
    fill.key = interaction(sw.family, outcome, sep = "_")
  )

## plot fraction of null hypo. rejections ----
ks.bar.aftcen.ex <- ks.bar.aftcen %>% 
  filter(filt == "ibdends", imap == "hapmap")
ggplot(
  ks.bar.aftcen,
  aes(x = param.code, y = frac, fill = fill.key)
) +
  geom_col(position = "stack", width = 0.8, color = "black", linewidth = 0.5) +
  facet_grid(imap ~ filt) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_manual(
    values = c(
      "hapibd_don't reject" = "lightblue1",
      "hapibd_reject"       = "darkblue",
      "rapid_don't reject"  = "mistyrose",
      "rapid_reject"        = "darkred"
    ),
    labels = c(
      "hapibd_don't reject" = "hap-IBD: don't reject",
      "hapibd_reject"       = "hap-IBD: reject",
      "rapid_don't reject"  = "RaPID: don't reject",
      "rapid_reject"        = "RaPID: reject"
    )
  ) +
  labs(
    title = glue(
      "K-S Test True vs Inferred {len.focal.pop} post-centromere IBD length distributions (unrelated pairs)"
    ),
    x = "Software / parameterization",
    y = "Fraction of replicates",
    fill = NULL
  ) +
  theme_bw(base_size = 32)
