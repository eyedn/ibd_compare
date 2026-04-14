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
auto.len <- sum(auto.lens$X2)

# calculate observed heterozygosity from observed homozygosity
mut.rate <- 1.3e-8
exp.het <- het.tbl %>%
  group_by(INDV) %>%
  summarise_at(.vars = c("N_SITES","O(HOM)"), .funs = sum) %>%
  ungroup() %>%
  mutate(`O(HET)`=((N_SITES-`O(HOM)`)/auto.len))
summary(exp.het$`O(HET)`)
shapiro.test(exp.het$`O(HET)`)

# plot mean and median heterozygotisties and Ne
mean.Ohet <- mean(exp.het$`O(HET)`)
med.Ohet <- median(exp.het$`O(HET)`)

mean.Ne <- mean.Ohet / (4 * mut.rate)
med.Ne <- med.Ohet / (4 * mut.rate)

ggplot(exp.het, aes(x = `O(HET)`)) +
  geom_histogram() +
  geom_vline(xintercept = mean.Ohet, color = "#D81B60") +
  geom_vline(xintercept = med.Ohet, color = "#1E88E5") +
  annotate(
    geom = "text",
    x = mean.Ohet,
    y = Inf,
    label = paste0("mean Ne = ", round(mean.Ne, 2)),
    color = "#D81B60",
    angle = 90,
    vjust = -0.5,
    hjust = 1.1
  ) +
  annotate(
    geom = "text",
    x = med.Ohet,
    y = Inf,
    label = paste0("median Ne = ", round(med.Ne, 2)),
    color = "#1E88E5",
    angle = 90,
    vjust = -0.5,
    hjust = 1.1
  ) +
  theme_bw()


# single ADMIX sim rep calculation ----
sim.het.tbl <- read_tsv(
  "~/scratch/stdpopsim_runs_chr18_admix2/vcfs/hetTest.het", col_names = TRUE, 
  show_col_types = FALSE
  )

chr18.len <- auto.lens$X2[18]
sim.mut.rate <- 2.36e-8

sim.exp.het <-  sim.het.tbl %>%
  group_by(INDV) %>%
  summarise_at(.vars = c("N_SITES","O(HOM)"), .funs = sum) %>%
  ungroup() %>%
  mutate(`O(HET)`=((N_SITES-`O(HOM)`)/chr18.len))
summary(sim.exp.het$`O(HET)`)
shapiro.test(sim.exp.het$`O(HET)`)

# plot mean and median heterozygotisties and Ne
sim.mean.Ohet <- mean(sim.exp.het$`O(HET)`)
sim.med.Ohet <- median(sim.exp.het$`O(HET)`)

sim.mean.Ne <- sim.mean.Ohet / (4 * sim.mut.rate)
sim.med.Ne <- sim.med.Ohet / (4 * sim.mut.rate)

ggplot(sim.exp.het, aes(x = `O(HET)`)) +
  geom_histogram() +
  geom_vline(xintercept = sim.mean.Ohet, color = "#D81B60") +
  geom_vline(xintercept = sim.med.Ohet, color = "#1E88E5") +
  annotate(
    geom = "text",
    x = sim.mean.Ohet,
    y = Inf,
    label = paste0("mean Ne = ", round(sim.mean.Ne, 2)),
    color = "#D81B60",
    angle = 90,
    vjust = -0.5,
    hjust = 1.1
  ) +
  annotate(
    geom = "text",
    x = sim.med.Ohet,
    y = Inf,
    label = paste0("median Ne = ", round(sim.med.Ne, 2)),
    color = "#1E88E5",
    angle = 90,
    vjust = -0.5,
    hjust = 1.1
  ) +
  theme_bw()
