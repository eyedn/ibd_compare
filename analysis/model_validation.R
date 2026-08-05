# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# model_validation.R
# ______________________________________________________________________________

library(tidyverse)

# helper functions ----
boot.mean.ci <- function(x, n.boot = 10000, conf = 0.95) {
  x <- x[is.finite(x)]
  stopifnot(length(x) > 0)

  boot.means <- replicate(n.boot, mean(sample(x, replace = TRUE)))
  alpha <- 1 - conf

  c(
    mean = mean(x),
    lo = quantile(boot.means, alpha / 2, names = FALSE),
    hi = quantile(boot.means, 1 - alpha / 2, names = FALSE)
  )
}

plot.ancestry.summary <- function(summary.tbl, value.col, ci, title, xlab) {
  legend.df <- tibble(
    xintercept = c(ci["mean"], ci["lo"], ci["hi"]),
    label = c(
      paste0("mean = ", round(ci["mean"], 4)),
      paste0("95% CI lower = ", round(ci["lo"], 4)),
      paste0("95% CI upper = ", round(ci["hi"], 4))
    ),
    type = c("Mean", "CI lower", "CI upper")
  )

  ggplot(summary.tbl, aes(x = {{ value.col }})) +
    geom_histogram() +
    geom_vline(
      data = legend.df,
      aes(xintercept = xintercept, color = type),
      linewidth = 1
    ) +
    scale_color_manual(
      values = c(
        "Mean" = "#D81B60",
        "CI lower" = "#1E88E5",
        "CI upper" = "#1E88E5"
      ),
      labels = legend.df$label,
      name = NULL
    ) +
    labs(
      x = xlab,
      y = "Count",
      title = title
    ) +
    theme_bw(base_size = 32) +
    theme(
      legend.position = "inside",
      legend.position.inside = c(0.02, 0.98),
      legend.justification = c(0, 1),
      legend.background = element_rect(fill = "white", color = "black"),
      legend.key = element_blank()
    )
}

# ADMIX het/Ne calculation from vcf het data ----
# read in GRCh39.p14 chromsome lengths
lens.file <- "~/references/GRCh38.p14/GRCh38.p14.chrlens.txt"
auto.lens <- read_tsv(lens.file, col_names = FALSE, show_col_types = FALSE) %>%
  filter(X1 %in% c(1:22)) 
tot.auto.len <- sum(auto.lens$X2)
chr18.len <- auto.lens$X2[18]

# read in heterozygosity results
sim.het.dir <- "~/scratch/msprime_runs2_chr18_admix_3T/ADMIX_het"
sim.het.tbl <- map_dfr(
  .x = 1:50,
  .f = ~ read_tsv(
    file.path(sim.het.dir, paste0(.x, ".het")),
    col_names = TRUE,
    show_col_types = FALSE
  ) %>%
    mutate(rep = .x)
)

# calculate observed heterozygosity from observed homozygosity
sim.mut.rate <- 2.36e-8
sim.het <- sim.het.tbl %>%
  mutate(`O(HET)`=((N_SITES-`O(HOM)`)/chr18.len))
summary(sim.het$`O(HET)`)
# shapiro.test(sim.het$`O(HET)`)

# plot mean and median heterozygotisties and Ne
sim.mean.Ohet <- mean(sim.het$`O(HET)`)
sim.med.Ohet <- median(sim.het$`O(HET)`)
sim.mean.Ne <- sim.mean.Ohet / (4 * sim.mut.rate)
sim.med.Ne <- sim.med.Ohet / (4 * sim.mut.rate)

sim.legend.df <- data.frame(
  xintercept = c(sim.mean.Ohet, sim.med.Ohet),
  label = c(
    paste0("mean Ne = ", round(sim.mean.Ne, 2)),
    paste0("median Ne = ", round(sim.med.Ne, 2))
  ),
  type = c("Mean", "Median")
)

ggplot(sim.het, aes(x = `O(HET)`)) +
  geom_histogram() +
  geom_vline(
    data = sim.legend.df,
    aes(xintercept = xintercept, color = type),
    linewidth = 1
  ) +
  scale_color_manual(
    values = c("Mean" = "#D81B60", "Median" = "#1E88E5"),
    labels = sim.legend.df$label,
    name = NULL
  ) +
  labs(
    x = "Observed heterozygosity",
    y = "Count",
    title = "Observed heterozygosity of ADMIX on chr18 (\"larger\" model)"
  ) +
  theme_bw(base_size = 32) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.02, 0.98), 
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_blank()
  )

# ADMIX ancestry proportions ----
anc.dir <- "~/scratch/msprime_runs_chr18_admix_3T/global_ancestry"

anc.files <- list.files(
  anc.dir,
  pattern = "_ADMIX\\.tsv$",
  full.names = TRUE
)
stopifnot(length(anc.files) > 0)

admix.anc <- map_dfr(
  anc.files,
  function(anc.file) {
    rep.id <- str_match(basename(anc.file), "_([0-9]+)_ADMIX\\.tsv$")[, 2]

    if (is.na(rep.id)) {
      stop("could not parse replicate ID from file: ", anc.file)
    }

    anc.tbl <- read_tsv(anc.file, col_names = TRUE, show_col_types = FALSE)

    tibble(
      rep = as.integer(rep.id),
      african.ancestry = anc.tbl[[2]]
    )
  }
)

admix.anc.summary <- admix.anc %>%
  group_by(rep) %>%
  summarise(
    mean.african.ancestry = mean(african.ancestry),
    sd.african.ancestry = sd(african.ancestry),
    .groups = "drop"
  )

mean.african.ancestry.ci <- boot.mean.ci(
  admix.anc.summary$mean.african.ancestry
)
sd.african.ancestry.ci <- boot.mean.ci(
  admix.anc.summary$sd.african.ancestry
)

plot.ancestry.summary(
  summary.tbl = admix.anc.summary,
  value.col = mean.african.ancestry,
  ci = mean.african.ancestry.ci,
  title = "Mean African ancestry of ADMIX (\"larger\" model)",
  xlab = "Replicate mean African ancestry"
)

plot.ancestry.summary(
  summary.tbl = admix.anc.summary,
  value.col = sd.african.ancestry,
  ci = sd.african.ancestry.ci,
  title = "SD African ancestry of ADMIX (\"larger\" model)",
  xlab = "Replicate SD African ancestry"
)

