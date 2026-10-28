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
library(glue)

# Resolve sources from either supported working directory.
ANALYSIS_DIR <- if (file.exists("analysis/R_utils/constants.R")) {
  "analysis"
} else "."
source(file.path(ANALYSIS_DIR, "R_utils/constants.R"), local = TRUE)
source(file.path(ANALYSIS_DIR, "R_utils/figure.core.R"), local = TRUE)
source(file.path(ANALYSIS_DIR, "R_utils/read.all.stats.R"), local = TRUE)
source(file.path(ANALYSIS_DIR, "R_utils/read.all.stats.per.rep.R"), local = TRUE)

# define constants ----
per.rep.path <- glue("{ibd.comp.dir}/true_vs_inference/all_stats.per_rep.tsv")

# read in replicate-level statistics ----
all.stats.per.rep <- read.all.stats.per.rep(
  per.rep.path,
  stats.meta.names,
  param.map
) %>%
  mutate(
    filt = factor(filt, levels = filt.levels),
    software.code = factor(software.code),
    pop = factor(pop)
  )


validate.figure.stats(all.stats.per.rep, ANALYSIS_PROFILES$code,
                      ANALYSIS_POPULATIONS, filt.levels, FOCAL_REPLICATE)

df.reltrue.all <- all.stats.per.rep %>%
  filter(imap == "hapmap") %>%
  pivot_longer(
    cols = matches("^ibd\\.score\\.rel\\.true(\\.|$)"),
    values_to = "ibd.score.rel.true"
  )

df.true.ibd.all <- all.stats.per.rep %>%
  filter(
    imap == "hapmap",
    rep == FOCAL_REPLICATE,
    filt == "ibdends",
    startsWith(as.character(software.code), "h")
  ) %>%
  pivot_longer(
    cols = matches("^ibd\\.score\\.err(\\.|$)"),
    values_to = "ibd.score.err"
  )
df.true.ibd.all <- df.true.ibd.all %>%
  group_by(software.code) %>%
  mutate(rel.score = safe.ratio(true.ibd.score,
    true.ibd.score[pop == "AFR"][1])) %>%
  ungroup()
df.err.all <- all.stats.per.rep %>%
  filter(
    imap == "hapmap",
    filt != "ibdendshot",
    # imap != "hapmap" | filt != "hot"
  ) %>%
  pivot_longer(
    cols = matches("^ibd\\.score\\.err(\\.|$)"),
    values_to = "ibd.score.err"
  )
df.err.afr.all <- all.stats.per.rep %>%
  # keep whatever imap/filt/software you want; don't force hapmap/ibdends if you want facets
  # filter(imap != "decode", filt != "none", imap != "hapmap" | filt != "hot") %>%
  filter(imap == "hapmap", filt == "ibdends") %>%
  group_by(imap, filt, software.code, rep) %>%   # add other ID cols here if needed (e.g., tmap)
  mutate(
    # pull the AFR baselines within each (imap,filt,software,rep) group
    true.afr = true.ibd.score[pop == "AFR"][1],
    inf.afr  = inf.ibd.score[pop == "AFR"][1],

    # normalize to AFR (relative-to-AFR scores)
    true.ibd.score.rel.afr = safe.ratio(true.ibd.score, true.afr),
    inf.ibd.score.rel.afr  = safe.ratio(inf.ibd.score, inf.afr),

    # error = (true relative) - (inferred relative)
    ibd.score.rel.afr.err = true.ibd.score.rel.afr - inf.ibd.score.rel.afr
  ) %>%
  ungroup()


ibd.score.plots <- list()
true.score.plots <- list()
for (pair.name in rownames(PROFILE_PAIRINGS)) {
  hapibd.code <- PROFILE_PAIRINGS[pair.name, "hapibd"]
  rapid.code <- PROFILE_PAIRINGS[pair.name, "rapid"]
  selected.codes <- c(hapibd.code, rapid.code)
  plots <- list()
  # analyze ibd score rel true performance ----
  df.reltrue <- df.reltrue.all %>% filter(software.code %in% selected.codes)
  plots[["true.box"]] <- ggplot(df.reltrue) +
    geom_boxplot(
      aes(x = pop, y = true.ibd.score),
      outliers = FALSE
    ) +
    facet_grid(filt ~ software, labeller = filt.labeller) +
    theme_bw(base_size = 32) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

  plots[["relative.true"]] <- ggplot(
    df.reltrue,
    aes(x = filt, y = ibd.score.rel.true, color = software.code)
  ) +
    geom_boxplot(
      aes(group = interaction(filt, software.code)), fill = "lightgrey"
    ) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_x_discrete(labels = filt.labels) +
    scale_color_manual(values = software.palette) +
    labs(
      x = "Filtering procedure",
      y = "Relative IBD score (inferred / true)",
      color = "Software",
      title = "Inferred vs True IBD score using ibd-ends"
    ) +
    coord_trans(y = "log10") +
    facet_wrap(~ pop) +
    theme_bw(base_size = 32)

  if (rapid.code == ANALYSIS_PROFILES$code[4]) {
  ## plot true ibd scores ----
  df.true.ibd <- df.true.ibd.all %>% filter(software.code == hapibd.code)
  plots[["true.score"]] <- ggplot(df.true.ibd, aes(x = pop, y = true.ibd.score, fill = pop)) +
    geom_col() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_manual(values = pop.palette) +
    labs(
      x = NULL,
      y = "True IBD score",
      title = glue("IBD score"),
      fill = NULL
    ) +
    theme_bw(base_size = 32) +
    theme(legend.position = "none")
  plots[["true.relative.afr"]] <- ggplot(df.true.ibd, aes(x = pop, y = rel.score, fill = pop)) +
    geom_col() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_manual(values = pop.palette) +
    labs(
      x = NULL,
      y = "True IBD score rel. to AFR",
      title = glue("IBD score rel. to AFR"),
      fill = NULL
    ) +
    theme_bw(base_size = 32) +
    theme(legend.position = "none")

    for (family in c("true.score", "true.relative.afr")) {
      plot <- plots[[family]]
      attr(plot, "profiles") <- hapibd.code
      true.score.plots[[family]][[paste(hapibd.code, "hapmap_ibdends",
        paste0("rep", FOCAL_REPLICATE), sep = "_")]] <- plot
      plots[[family]] <- NULL
    }
  }
  # analyze ibd score error performance ----
  df.err <- df.err.all %>% filter(software.code %in% selected.codes)
  plots[["score.error"]] <- ggplot(df.err, aes(x = pop, y = ibd.score.err, fill = software.code)) +
    geom_boxplot() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_manual(values = software.palette) +
    labs(
      x = NULL,
      y = "IBD score error",
      title = glue("IBD score error"),
      fill = NULL
    ) +
    facet_grid(filt ~ software, labeller = filt.labeller) +
    theme_bw(base_size = 32)

  # plot rel ibd score error to AFR ----
  df.err <- df.err.afr.all %>% filter(software.code %in% selected.codes)
  plots[["relative.afr.error"]] <- ggplot(df.err, aes(x = pop, y = ibd.score.rel.afr.err, fill = software.code)) +
    geom_boxplot(outlier.shape = NA) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_manual(values = software.palette) +
    labs(
      x = NULL,
      y = "Relative IBD score error",
      title = "IBD score error (relative to AFR)" ,
      fill = NULL
    ) +
    theme_bw(base_size = 32)

  ibd.score.plots[[pair.name]] <- plots
}
ibd.score.plots <- c(organize.figure.plots(ibd.score.plots), true.score.plots)
invisible(ibd.score.plots)
