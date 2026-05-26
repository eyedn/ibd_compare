# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# model.validation.R
# ______________________________________________________________________________

library(tidyverse)

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
shapiro.test(sim.het$`O(HET)`)

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

# ADMIX ancestry proportions ----
#TODO: read in all replicate global anc. prop. files
anc.dir <- "~/scratch/msprime_runs_chr18_admix_3T/global_ancestry"

#TODO: get CI for mean african ancestry; african ancestry stored in X2
#TODO: plot above distribution with CI and mean mean

#TODO: get CI for sd african ancestry; african ancestry stored in X2
#TODO: plot above distribution with CI and mean sd
