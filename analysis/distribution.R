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
library(glue)

# Resolve sources from either supported working directory.
ANALYSIS_DIR <- if (file.exists("analysis/R_utils/constants.R")) {
  "analysis"
} else "."
source(file.path(ANALYSIS_DIR, "R_utils/constants.R"), local = TRUE)
source(file.path(ANALYSIS_DIR, "R_utils/figure.core.R"), local = TRUE)
source(file.path(ANALYSIS_DIR, "R_utils/read.figure.inputs.R"), local = TRUE)
source(file.path(ANALYSIS_DIR, "R_utils/read.pop.calls.R"), local = TRUE)
source(file.path(ANALYSIS_DIR, "R_utils/read.rep.calls.R"), local = TRUE)
source(file.path(ANALYSIS_DIR, "R_utils/run.ks.R"), local = TRUE)

cores <- WORKERS
unrels.tbl <- read.population.files(POP_INFO_DIR, CHROMOSOME)
len.focal.pop <- "EUR"
len.focal.rep <- FOCAL_REPLICATE
imap.all.calls <- read.pop.calls(
  len.focal.pop, dir.paths.map, unrels.tbl, cores
)
rep1.calls <- imap.all.calls %>% filter(rep == len.focal.rep)
# analyze length distribution using ks testing ----
## ks testing procedure ----
ks.res <- ks.figure.summary(imap.all.calls, attr(imap.all.calls, "datasets"))

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

# analyze length distribution post-centromere using ks testing ----
## ks testing procedure ----
ks.res.aftcen <- ks.figure.summary(
  imap.all.calls, attr(imap.all.calls, "datasets"), after.centromere = TRUE
)

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

ks.all <- ks.bar
ks.aftcen.all <- ks.bar.aftcen
distribution.plots <- list()
distribution.single <- list()
for (pair.name in rownames(PROFILE_PAIRINGS)) {
  hapibd.code <- PROFILE_PAIRINGS[pair.name, "hapibd"]
  rapid.code <- PROFILE_PAIRINGS[pair.name, "rapid"]
  selected.codes <- c(hapibd.code, rapid.code)
  ks.bar <- ks.all %>% filter(param.code %in% selected.codes)
  ks.bar.aftcen <- ks.aftcen.all %>% filter(param.code %in% selected.codes)
  plots <- list()
  ## plot fraction of null hypo. rejections ----
  ks.bar.ex <- ks.bar %>% filter(filt == "ibdends", imap == "hapmap")
  plots[["ks.fraction"]] <- ggplot(
    ks.bar,
    aes(x = param.code, y = frac, fill = fill.key)
  ) +
    geom_col(position = "stack", width = 0.8, color = "black", linewidth = 0.5) +
    facet_grid(imap ~ filt, labeller = filt.labeller) +
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

  if (hapibd.code == ANALYSIS_PROFILES$code[1]) {
  # plot ibd distributions of rep 1 ----
  ## plot ibd length distribution ----
  sub.call.y.max <- 500 # specific to each pop
  sub.calls <- rep1.calls %>%
    filter(
      param.code == "true" |
        (param.code == rapid.code & filt == "ibdends" & imap == "hapmap") |
        (param.code == rapid.code & filt == "ibdends" & imap == "1Mb1cM"),
      pop == len.focal.pop,
      is.unrel == TRUE
    ) %>%
    # hard coded facet wrap order
    mutate(
      filt.label = filt.labels[as.character(filt)],
      facet = interaction(param.code, imap, filt.label, sep = ", "),
      facet = fct_relevel(
        facet,
        "true, true, N",
        paste0(rapid.code, ", 1Mb1cM, I"),
        paste0(rapid.code, ", hapmap, I")
      )
    )

  plots[["length"]] <- ggplot(sub.calls) +
    geom_histogram(
      aes(x = cM.len, fill = param.code), color = "black",  linewidth = 0.5
    ) +
    facet_wrap(~ facet, drop = TRUE) +
    coord_cartesian(ylim = c(0, sub.call.y.max)) +
    scale_fill_manual(
      values = software.palette,
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
      fill.plot = factor(fill.plot, levels = c("true", hapibd.code, rapid.code, "cen"))
    )

  plots[["centromere"]] <- ggplot(sub.calls.cen) +
    geom_histogram(
      aes(x = cM.len, fill = fill.plot),
      position = "stack",
      color = "black",
      linewidth = 0.5
    ) +
    facet_wrap(~ facet, drop = TRUE) +
    coord_cartesian(ylim = c(0, sub.call.y.max)) +
    scale_fill_manual(
      values = c(software.palette, "true" = "lightgrey", "cen" = "black")
    ) +
    labs(
      title = glue(
        "{len.focal.pop} IBD segments overlapping centromere (black bars)"
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
      fill.plot = factor(fill.plot, levels = c("true", hapibd.code, rapid.code, "bfcen"))
    )

  plots[["before.centromere"]] <- ggplot(sub.calls.bfcen) +
    geom_histogram(
      aes(x = cM.len, fill = fill.plot),
      position = "stack",
      color = "black",
      linewidth = 0.5
    ) +
    facet_wrap(~ facet, drop = TRUE) +
    coord_cartesian(ylim = c(0, sub.call.y.max)) +
    scale_fill_manual(
      values = c(software.palette, "true" = "lightgrey", "bfcen" = "black")
    ) +
    labs(
      title = glue(
        "{len.focal.pop} IBD segments before centromere (black bars)"
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
      fill.plot = factor(fill.plot, levels = c("true", hapibd.code, rapid.code, "aftcen"))
    )

  plots[["after.centromere"]] <- ggplot(sub.calls.aftcen) +
    geom_histogram(
      aes(x = cM.len, fill = fill.plot),
      position = "stack",
      color = "black",
      linewidth = 0.5
    ) +
    facet_wrap(~ facet, drop = TRUE) +
    coord_cartesian(ylim = c(0, sub.call.y.max)) +
    scale_fill_manual(
      values = c(software.palette, "true" = "lightgrey", "aftcen" = "black")
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

    for (family in c("length", "centromere", "before.centromere",
                     "after.centromere")) {
      plot <- plots[[family]]
      attr(plot, "profiles") <- rapid.code
      distribution.single[[family]][[paste(rapid.code, len.focal.pop,
        paste0("rep", FOCAL_REPLICATE), sep = "_")]] <- plot
      plots[[family]] <- NULL
    }
  }
  ## plot fraction of null hypo. rejections ----
  ks.bar.aftcen.ex <- ks.bar.aftcen %>%
    filter(filt == "ibdends", imap == "hapmap")
  plots[["ks.post.centromere"]] <- ggplot(
    ks.bar.aftcen,
    aes(x = param.code, y = frac, fill = fill.key)
  ) +
    geom_col(position = "stack", width = 0.8, color = "black", linewidth = 0.5) +
    facet_grid(imap ~ filt, labeller = filt.labeller) +
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

  # figure for PEQG ----
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
  ks.df <- ks.bar %>%
    filter(imap != "decode", param.code %in% selected.codes) %>%
    mutate(perc = frac * 100)
  ks.plot <- ggplot(
    ks.df,
    aes(x = param.code, y = perc, fill = fill.key)
  ) +
    geom_col(position = "stack", width = 0.8, color = "black", linewidth = 1) +
    facet_grid(
      imap ~ filt,
      labeller = labeller(imap = map.labels.tmp, filt = filt.labels.tmp)
    ) +
    scale_y_continuous(
      breaks = c(0, 25, 50, 75, 100)
    ) +
    scale_fill_manual(
      values = c(
        "hapibd_don't reject" = "#EBC9DE",
        "hapibd_reject"       = "#C75DAA",
        "rapid_don't reject"  = "#B0DDDE",
        "rapid_reject"        = "#009B9F"
      ),
      labels = c(
        "hapibd_don't reject" = "hap-IBD: don't reject",
        "hapibd_reject"       = "hap-IBD: reject",
        "rapid_don't reject"  = "RaPID: don't reject",
        "rapid_reject"        = "RaPID: reject"
      )
    ) +
    labs(
      title =  NULL,
      x = NULL,
      y = "% replicates",
      fill = NULL
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


  # save precision figure for PEQG

  plots$ks.presentation <- ks.plot
  distribution.plots[[pair.name]] <- plots
}
distribution.plots <- c(organize.figure.plots(distribution.plots),
                        distribution.single)
invisible(distribution.plots)
