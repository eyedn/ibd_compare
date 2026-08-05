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
library(furrr)
library(future)
library(progressr)
library(glue)
library(stats)
library(emmeans)

# source function ----
source("R_utils/constants.R")
source("R_utils/read.plink2.king.R")
source("R_utils/kinship.mat2tbl.R")
source("R_utils/read.pop.ibd.kinship.R")
source("R_utils/read.rep.ibd.kinship.R")

# set up for parallel processing ----
cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
handlers("cli")

# read in KING kinships ----
king.dir <- glue("{sim.dir}/pop_info/plink2_king_out")
king.tbl <- tibble(
  file = list.files(
    king.dir, full.names = TRUE,
    pattern = "^HapMapII_GRCh38_\\d+_[^_]+\\.kin0$"
  )
) %>%
  mutate(
    rep = as.integer(str_match(basename(file), "^HapMapII_GRCh38_(\\d+)_")[, 2]),
    pop = str_match(basename(file), "^HapMapII_GRCh38_\\d+_([^_]+)\\.kin0$")[, 2],
    data = map(file, read.plink2.king)
  ) %>%
  select(rep, pop, data) %>%
  unnest(data)
king.tbl <- king.tbl %>%
  mutate(
    sam1.min = pmin(sam1, sam2),
    sam2.max = pmax(sam1, sam2),
    sam.pair = paste0(sam1.min, "-", sam2.max)
  ) %>%
  select(-sam1.min, -sam2.max)

# read in kinship results for a specific pop ----
kin.focal.pop <- "ADX"
true.tbl <- read.pop.ibd.kinship(
  kin.focal.pop, tvi.dir.paths.map, cores, which = "true"
  ) %>%
  select(-filt, -software.code)
inf.tbl <- read.pop.ibd.kinship(
  kin.focal.pop, tvi.dir.paths.map, cores, which = "inf"
  )
inf.tbl <- list(
    true.tbl %>% rename(kin.true = kinship) %>% select(-sam1, -sam2),
    inf.tbl  %>% rename(kin.inf  = kinship) %>% select(-sam1, -sam2)
  ) %>%
  reduce(full_join, by = c("rep", "pop", "sam.pair"))
summary(inf.tbl$kin.true)
tmp <- inf.tbl %>% filter(software.code == "rf1", filt == "none")
cor.test(tmp$kin.true, tmp$kin.inf)

# read in kinships result for rep 1 ----
kin.focal.rep <- 1
true.rep1.tbl <- read.rep.ibd.kinship(
    kin.focal.rep, tvi.dir.paths.map, cores, which = "true"
  ) %>%
  select(-filt, -software.code)
inf.rep1.tbl <- read.rep.ibd.kinship(
    kin.focal.rep, tvi.dir.paths.map, cores, which = "inf"
  )
inf.rep1.tbl <- list(
  true.rep1.tbl %>% rename(kin.true = kinship) %>% select(-sam1, -sam2),
  inf.rep1.tbl  %>% rename(kin.inf  = kinship) %>% select(-sam1, -sam2)
) %>%
  reduce(full_join, by = c("rep", "pop", "sam.pair"))

# analyze kindship performance ----
# ---- set these two dirs correctly ----
hap.dir <- "~/scratch/ibd_compare_chr18_admix_3T/true_vs_inference/rapid_sim_tmapHapMapII_GRCh38_imapHapMapII_GRCh38_d2_w100_r10_s2/ibd_kinship_mats"

# ---- read mats ----
true.hap <- readRDS(file.path(hap.dir, "true_ADX_30.ibdends.rds"))
inf.hap  <- readRDS(file.path(hap.dir, "inf_ADX_30.ibdends.rds"))

true.tbl <- kinship.mat2tbl(true.hap, 30, "ADX")
inf.tbl <- kinship.mat2tbl(inf.hap, 30, "ADX")
king.sub.tbl <- king.tbl %>% filter(rep == "30", pop == "ADX")

joined.tbl <- list(
  true.tbl %>% rename(kin.true = kinship) %>% select(-sam1, -sam2),
  inf.tbl  %>% rename(kin.inf  = kinship) %>% select(-sam1, -sam2),
  king.sub.tbl %>% rename(kin.king = kinship) %>% select(-sam1, -sam2)
) %>%
  reduce(full_join, by = c("rep", "pop", "sam.pair")) %>%
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
    by = c("rep", "pop", "sam.pair")
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
    r = cor(kin.x, kin.y, use = "complete.obs"),
    .groups = "drop"
  ) %>%
  mutate(label = sprintf("r = %.3f", r))

ggplot(plot.pairs, aes(x = kin.x, y = kin.y)) +
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

# results for PEQG ----
# read in data
true.ADMIX.tbl <- read.pop.ibd.kinship(
  "ADX", tvi.dir.paths.map, cores, which = "true"
) %>%
  select(-filt, -software.code)
inf.ADMIX.tbl <- read.pop.ibd.kinship(
  "ADX", tvi.dir.paths.map, cores, which = "inf"
)
inf.ADMIX.tbl <- list(
  true.ADMIX.tbl %>% rename(kin.true = kinship) %>% select(-sam1, -sam2),
  inf.ADMIX.tbl  %>% rename(kin.inf  = kinship) %>% select(-sam1, -sam2)
) %>%
  reduce(full_join, by = c("rep", "pop", "sam.pair"))

true.AFR.tbl <- read.pop.ibd.kinship(
  "AFR", tvi.dir.paths.map, cores, which = "true"
) %>%
  select(-filt, -software.code)
inf.AFR.tbl <- read.pop.ibd.kinship(
  "AFR", tvi.dir.paths.map, cores, which = "inf"
)
inf.AFR.tbl <- list(
  true.AFR.tbl %>% rename(kin.true = kinship) %>% select(-sam1, -sam2),
  inf.AFR.tbl  %>% rename(kin.inf  = kinship) %>% select(-sam1, -sam2)
) %>%
  reduce(full_join, by = c("rep", "pop", "sam.pair"))

true.EUR.tbl <- read.pop.ibd.kinship(
  "EUR", tvi.dir.paths.map, cores, which = "true"
) %>%
  select(-filt, -software.code)
inf.EUR.tbl <- read.pop.ibd.kinship(
  "EUR", tvi.dir.paths.map, cores, which = "inf"
)
inf.EUR.tbl <- list(
  true.EUR.tbl %>% rename(kin.true = kinship) %>% select(-sam1, -sam2),
  inf.EUR.tbl  %>% rename(kin.inf  = kinship) %>% select(-sam1, -sam2)
) %>%
  reduce(full_join, by = c("rep", "pop", "sam.pair"))

inf.tbl <- bind_rows(inf.ADMIX.tbl, inf.AFR.tbl, inf.EUR.tbl)
inf.tbl <- inf.tbl %>%
  mutate(
    imap = ifelse(substr(software.code, 2, 2) == "f", "flat", "HapMapII"),
    software.code = paste0(
      substr(software.code, 1, 1), 
      substr(
        software.code, nchar(software.code), nchar(software.code)
        )
      )
    )
inf.tbl <- inf.tbl %>% filter(kin.true > 0)

# data formatting
kin.summary <- inf.tbl %>%
  group_by(pop, filt, imap, software.code) %>%
  summarise(
    n = n(),
    pearson_r = cor(
      kin.true,
      kin.inf,
      method = "pearson",
      use = "complete.obs"
    ),
    spearman_rho = cor(
      kin.true,
      kin.inf,
      method = "spearman",
      use = "complete.obs"
    ),
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
    values = c(
      "h3" = "#C75DAA",
      "r1" = "#009B9F"
    ),
    labels = c(
      "h3"   = "hap-IBD",
      "r1"   = "RaPID"
    )
  ) +
  scale_fill_manual(
    name = NULL,
    values = c(
      "h3"   = "#EBC9DE",
      "r1"   = "#B0DDDE"
    ),
    labels = c(
      "h3"   = "hap-IBD",
      "r1"   = "RaPID"
    )
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

if (interactive()) {
  print(kin.plot)
}

pdf(
  "~/proj/ibd_compare_results/kin_flat.pdf",
  width = 10,
  height = 5,
  family = "Helvetica"
)

print(kin.plot)

dev.off()
