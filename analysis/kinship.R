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
kin.focal.pop <- "ADMIX"
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
hap.dir <- "~/scratch/ibd_compare_chr18_admix/true_vs_inference/rapid_sim_tmapHapMapII_GRCh38_imapHapMapII_GRCh38_d2_w100_r10_s2/ibd_kinship_mats"

# ---- read mats ----
true.hap <- readRDS(file.path(hap.dir, "true_ADMIX_30.ibdends.rds"))
inf.hap  <- readRDS(file.path(hap.dir, "inf_ADMIX_30.ibdends.rds"))

true.tbl <- kinship.mat2tbl(true.hap, 30, "ADMIX")
inf.tbl <- kinship.mat2tbl(inf.hap, 30, "ADMIX")
king.sub.tbl <- king.tbl %>% filter(rep == "30", pop == "ADMIX")

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

