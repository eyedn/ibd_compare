# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# get.ASWNe.R
# ______________________________________________________________________________


library(tidyverse)

# ASW calculation ----
# read in heterozygosity results
het.dir <- "~/proj/1000GenomeNYGC_hg38_karatas/ASW_het"
het.tbl <- map_dfr(
  .x = 1:22,
  .f = ~ read_tsv(
    file.path(het.dir, paste0(.x, ".het")),
    col_names = TRUE,
    show_col_types = FALSE
  ) %>%
    mutate(chr = .x)
)

# read in GRCh39.p14 chromsome lengths
lens.file <- "~/references/GRCh38.p14/GRCh38.p14.chrlens.txt"
auto.lens <- read_tsv(lens.file, col_names = FALSE, show_col_types = FALSE) %>%
  filter(X1 %in% c(1:22)) 
tot.auto.len <- sum(auto.lens$X2)
chr18.len <- auto.lens$X2[18]

# read in GRCh39.p14 gencode intergenic region lengths
intergenic.file <- "~/references/GRCh38.p14/gencodeV49.canonicalIntergenic.bed"
intergentic.regions <- read_tsv(
  intergenic.file, col_names = FALSE, show_col_types = FALSE
  )
intergentic.regions <- intergentic.regions %>%
  mutate(len = X3 - X2)
mut.rate <- 1.3e-8 

## all ASW chr ----
# calculate observed heterozygosity from observed homozygosity
tot.inter.len <- sum(intergentic.regions$len)
tot.het <- het.tbl %>%
  group_by(INDV) %>%
  summarise_at(.vars = c("N_SITES","O(HOM)"), .funs = sum) %>%
  ungroup() %>%
  mutate(`O(HET)`=((N_SITES-`O(HOM)`)/tot.auto.len))
summary(tot.het$`O(HET)`)
shapiro.test(tot.het$`O(HET)`)

# plot mean and median heterozygotisties and Ne
tot.mean.Ohet <- mean(tot.het$`O(HET)`)
tot.med.Ohet <- median(tot.het$`O(HET)`)
tot.mean.Ne <- tot.mean.Ohet / (4 * mut.rate)
tot.med.Ne <- tot.med.Ohet / (4 * mut.rate)

tot.legend.df <- data.frame(
  xintercept = c(tot.mean.Ohet, tot.med.Ohet),
  label = c(
    paste0("mean Ne = ", round(tot.mean.Ne, 2)),
    paste0("median Ne = ", round(tot.med.Ne, 2))
  ),
  type = c("Mean", "Median")
)

ggplot(tot.het, aes(x = `O(HET)`)) +
  geom_histogram() +
  geom_vline(
    data = tot.legend.df,
    aes(xintercept = xintercept, color = type),
    linewidth = 1
  ) +
  scale_color_manual(
    values = c("Mean" = "#D81B60", "Median" = "#1E88E5"),
    labels = tot.legend.df$label,
    name = NULL
  ) +
  labs(
    x = "Observed heterozygosity",
    y = "Count",
    title = "Observed heterozygosity of AFR"
  ) +
  theme_bw(base_size = 32) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.02, 0.98), 
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_blank()
  )

## just ASW chr18 ----
# calculate observed heterozygosity from observed homozygosity
chr18.inter.len <- sum(
  intergentic.regions[intergentic.regions$X1=="chr18", "len"]
  )
chr18.het <- het.tbl %>%
  group_by(INDV) %>%
  summarise_at(.vars = c("N_SITES","O(HOM)"), .funs = sum) %>%
  ungroup() %>%
  mutate(`O(HET)`=((N_SITES-`O(HOM)`)/chr18.inter.len))
summary(chr18.het$`O(HET)`)
shapiro.test(chr18.het$`O(HET)`)

# plot mean and median heterozygotisties and Ne
chr18.mean.Ohet <- mean(chr18.het$`O(HET)`)
chr18.med.Ohet <- median(chr18.het$`O(HET)`)
chr18.mean.Ne <- chr18.mean.Ohet / (4 * mut.rate)
chr18.med.Ne <- chr18.med.Ohet / (4 * mut.rate)

chr18.legend.df <- data.frame(
  xintercept = c(chr18.mean.Ohet, chr18.med.Ohet),
  label = c(
    paste0("mean Ne = ", round(chr18.mean.Ne, 2)),
    paste0("median Ne = ", round(chr18.med.Ne, 2))
  ),
  type = c("Mean", "Median")
)

ggplot(chr18.het, aes(x = `O(HET)`)) +
  geom_histogram() +
  geom_vline(
    data = chr18.legend.df,
    aes(xintercept = xintercept, color = type),
    linewidth = 1
  ) +
  scale_color_manual(
    values = c("Mean" = "#D81B60", "Median" = "#1E88E5"),
    labels = chr18.legend.df$label,
    name = NULL
  ) +
  labs(
    x = "Observed heterozygosity",
    y = "Count",
    title = "Observed heterozygosity of AFR on chr18"
  ) +
  theme_bw(base_size = 32) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.02, 0.98), 
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_blank()
  )

# ADMIX calculation ----
# read in heterozygosity results
sim.het.dir <- "~/scratch/stdpopsim_runs_chr18_admix2/ADMIX_het"
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
shapiro.test(sim.het$`O(HET)`)

# plot mean and median heterozygotisties and Ne
sim.mean.Ohet <- mean(sim.het$`O(HET)`)
sim.med.Ohet <- median(sim.het$`O(HET)`)
sim.mean.Ne <- sim.mean.Ohet / (4 * mut.rate)
sim.med.Ne <- sim.med.Ohet / (4 * mut.rate)

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
    title = "Observed heterozygosity of ADMIX on chr18"
  ) +
  theme_bw(base_size = 32) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.02, 0.98), 
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_blank()
  )

