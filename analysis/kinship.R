# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# kinship.R
# ______________________________________________________________________________


library(tidyverse)
library(glue)

ANALYSIS_DIR <- if (file.exists("analysis/R_utils/constants.R")) {
  "analysis"
} else "."
source(file.path(ANALYSIS_DIR, "R_utils/constants.R"), local = TRUE)
source(file.path(ANALYSIS_DIR, "R_utils/figure.core.R"), local = TRUE)
source(file.path(ANALYSIS_DIR, "R_utils/read.plink2.king.R"), local = TRUE)
source(file.path(ANALYSIS_DIR, "R_utils/kinship.mat2tbl.R"), local = TRUE)
source(file.path(ANALYSIS_DIR, "R_utils/read.figure.inputs.R"), local = TRUE)
source(file.path(ANALYSIS_DIR, "R_utils/read.pop.ibd.kinship.R"), local = TRUE)
source(file.path(ANALYSIS_DIR, "R_utils/read.rep.ibd.kinship.R"), local = TRUE)

# Read each dataset once, retaining all maps and parameter strings.
cores <- WORKERS
king.tbl <- read.population.files(POP_INFO_DIR, CHROMOSOME, king = TRUE) %>%
  mutate(sam.pair = paste0(pmin(sam1, sam2), "-", pmax(sam1, sam2)))
true.tbl <- read.figure.kinship(tvi.dir.paths.map, cores, which = "true") %>%
  select(rep, pop, sam.pair, kin.true = kinship)
inf.tbl <- read.figure.kinship(tvi.dir.paths.map, cores, which = "inf") %>%
  rename(kin.inf = kinship) %>%
  left_join(true.tbl, by = c("rep", "pop", "sam.pair"))
joined.all <- inf.tbl

# The historical scatter family shows RaPID only, once per profile.
plots <- list()
for (rapid.code in ANALYSIS_PROFILES$code[4:6]) {
  joined.tbl <- joined.all %>%
    filter(pop == "ADX", rep == FOCAL_REPLICATE, imap == "hapmap",
           filt == "ibdends", param.code == rapid.code) %>%
    select(rep, pop, sam.pair, kin.true, kin.inf) %>%
    full_join(king.tbl %>%
      filter(pop == "ADX", rep == FOCAL_REPLICATE) %>%
      select(rep, pop, sam.pair, kin.king = kinship),
      by = c("rep", "pop", "sam.pair")) %>%
    mutate(kin.king.trunc = pmax(kin.king, 0))
# long table
kin.long <- joined.tbl %>%
  pivot_longer(
    cols = c(kin.true, kin.inf, kin.king, kin.king.trunc),
    names_to = "method",
    values_to = "kin"
  ) %>%
  mutate(
    method = recode(
      method,
      kin.true = "True IBD",
      kin.inf = "Inf. IBD",
      kin.king = "KING",
      kin.king.trunc = "KING trunc."
    )
  )

# all method pairs for each (rep,pop,sam.pair)
plot.pairs <- kin.long %>%
  rename(method.x = method, kin.x = kin) %>%
  inner_join(
    kin.long %>% rename(method.y = method, kin.y = kin),
    by = c("rep", "pop", "sam.pair"), relationship = "many-to-many"
  ) %>%
  filter(method.x != method.y) %>%
  # keep only the 6 panels you want
  filter(method.x %in% c("True IBD", "Inf. IBD")) %>%
  filter(!(method.x == "True IBD" & method.y == "Inf. IBD")) %>%
  filter(method.x != method.y)

# correlation labels per facet
cor.lbl <- plot.pairs %>%
  group_by(method.x, method.y) %>%
  summarise(
    r = safe.cor(kin.x, kin.y),
    .groups = "drop"
  ) %>%
  mutate(label = sprintf("r = %.3f", r))

plots[[paste(rapid.code, "ADX", paste0("rep", FOCAL_REPLICATE),
                "hapmap_ibdends", sep = "_")]] <- ggplot(plot.pairs, aes(x = kin.x, y = kin.y)) +
  geom_point(alpha = 0.15, size = 0.6) +
  geom_abline(linetype = "dashed") +
  geom_smooth(method = "lm", se = FALSE) +
  geom_text(
    data = cor.lbl,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = -0.05, vjust = 1.05,
    size = 6
  ) +
  facet_wrap( ~ method.y + method.x) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 32)

}
for (plot.name in names(plots)) {
  attr(plots[[plot.name]], "profiles") <- sub("_.*", "", plot.name)
}
kinship.plots <- list(scatter = plots)

# Keep the existing positive-truth restriction and flat-map summary figure.
inf.tbl <- joined.all %>%
  mutate(software.code = param.code,
         imap = recode(imap, hapmap = "HapMapII", `1Mb1cM` = "flat",
                       decode = "deCODE")) %>%
  filter(kin.true > 0)
# data formatting
kin.summary <- inf.tbl %>%
  group_by(pop, filt, imap, software.code) %>%
  summarise(
    n = n(),
    pearson_r = safe.cor(kin.true, kin.inf, method = "pearson"),
    spearman_rho = safe.cor(kin.true, kin.inf, method = "spearman"),
    rmse = sqrt(mean((kin.inf - kin.true)^2)),
    bias = mean(kin.inf - kin.true),
    .groups = "drop"
  )

# cor.tests <- inf.tbl %>%
#   group_by(pop, filt, imap, software.code) %>%
#   summarise(
#     test = list(cor.test(kin.true, kin.inf)),
#     .groups = "drop"
#   ) %>%
#   mutate(tidy = map(test, broom::tidy)) %>%
#   tidyr::unnest(tidy)


summary.all <- kin.summary
pair.plots <- list()
line.sizes <- 2
filt.labels.tmp <- c(none = "None", hot = "Z", ibdends = "I",
                     hotibdends = "Z then I", ibdendshot = "I then Z")
map.labels.tmp <- c(`1Mb1cM` = "flat", hapmap = "HapMapII")
for (pair.name in rownames(PROFILE_PAIRINGS)) {
  selected.codes <- unlist(PROFILE_PAIRINGS[pair.name, ], use.names = FALSE)
  kin.summary <- summary.all %>% filter(software.code %in% selected.codes)
  # kinships plot
  kin.summary <- kin.summary %>%
    mutate(
      imap = factor(imap, c("HapMapII", "flat")),
      filt = factor(filt, c("none", "hot", "ibdends", "hotibdends", "ibdendshot")),
      pop = factor(pop, c("ADX", "AFR", "EUR"))
      ) %>%
    filter(imap == "flat")

  kin.plot <- ggplot(
    kin.summary,
    aes(x = interaction(software.code, pop), y = spearman_rho)
  ) +
    geom_point(
      aes(
        color = software.code,
        fill = software.code,
        shape = pop
      ),
      size = 5, stroke = 2
    ) +
    facet_grid(
      . ~ filt,
      labeller = labeller(
        imap = map.labels.tmp,
        filt = filt.labels.tmp
      )
    ) +
    scale_shape_manual(
      name = NULL,
      values = c(
        "ADX" = 21,
        "AFR" = 24,
        "EUR" = 25
      ),
      labels = c(
        "ADX" = "ADX",
        "AFR" = "AFR",
        "EUR" = "EUR"
      )
    ) +
    scale_color_manual(
      name = NULL,
      values = setNames(c("#C75DAA", "#009B9F"), selected.codes),
      labels = setNames(c("hap-IBD", "RaPID"), selected.codes)
    ) +
    scale_fill_manual(
      name = NULL,
      values = setNames(c("#EBC9DE", "#B0DDDE"), selected.codes),
      labels = setNames(c("hap-IBD", "RaPID"), selected.codes)
    ) +
    labs(y = "flat Spearman's") +
    scale_y_continuous(limits = c(0.15, 1.05)) +
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
      shape = guide_legend(
        override.aes = list(
          fill = "white",
          color = "black"
        )
      ),
      fill = guide_legend(
        nrow = 1,
        byrow = TRUE,
        override.aes = list(
          shape = 21,
          color = c("#C75DAA", "#009B9F"),
          fill = c("#EBC9DE", "#B0DDDE"),
          size = 5,
          stroke = 2
        )
      )
    )

  pair.plots[[pair.name]] <- list(spearman.flat = kin.plot)
}
kinship.plots <- c(kinship.plots, organize.figure.plots(pair.plots))
invisible(kinship.plots)
