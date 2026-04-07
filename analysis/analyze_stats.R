# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology 
# Mooney Lab
# ___
# analyze_stats.R
# ______________________________________________________________________________


library(tidyverse)
library(furrr)
library(future)
library(progressr)
library(glue)
library(stats)
library(emmeans)

# # source function ----
# source("analyze_stats/rank.params.R")
# source("analyze_stats/run.ks.R")
# source("analyze_stats/compute.kinship.coeffs.R")
# source("analyze_stats/read.plink2.king.R")
# source("analyze_stats/constants.R")
# 
# # set up for parallel processing ----
# cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
# handlers("cli")
# 
# # read in statistics aggregated across replicates ----
# all.stats <- readr::read_tsv(
#   glue("{ibd.comp.dir}/true_vs_inference/all_stats.tsv"), 
#   col_names = TRUE, show_col_types = FALSE
#   )
# all.stats <- all.stats %>%
#   mutate(
#     software.code = param.map[software.with.params]
#   )
# 
# # read in statistics per replicates ----
# base.metrics <- c(
#   "tp","fn","fp","tn","sensitivity","specificity","ppv","npv","fpr","accuracy"
#   )
# suffixes <- c("", ".2", ".2.5", ".3", ".3.5", ".4", ".4.5", ".5")
# metric.names <- unlist(
#   lapply(suffixes, function(suf) paste0(base.metrics, suf)), use.names = FALSE
#   )
# meta.names <- c(
#   "software", "software.with.params", "tmap", "imap", "filt", "pop", "rep",
#   "n.true.unrels", "n.inf.unrels", "true.ibd.score", "inf.ibd.score",
#   "ibd.score.rel.true", "ibd.score.err"
# )
# all.stats.per.rep <- readr::read_tsv(
#   glue("{ibd.comp.dir}/true_vs_inference/all_stats.per_rep.tsv"), 
#   col_names = FALSE, show_col_types = FALSE
# )
# stopifnot(ncol(all.stats.per.rep) >= length(meta.names))
# n_metrics <- ncol(all.stats.per.rep) - length(meta.names)
# colnames(all.stats.per.rep) <- c(meta.names, metric.names)
# all.stats.per.rep <- all.stats.per.rep %>%
#   mutate(
#     software.code = param.map[software.with.params]
#   )

# read genetic maps ----
hapmap <- read_table(
  glue("{ref.dir}/HapMapII_GRCh38/plink.chr18.GRCh38.UCSC.format.map"), 
  col_names = FALSE, show_col_types = FALSE
  ) %>% rename(chr = X1, snp = X2, cM = X3, bp = X4)
decode <- read_tsv(
  glue("{ref.dir}/deCODE/addendum/sex-averaged_gmaps/sex_av_chr18_hg38.gmap"), 
  col_names = FALSE, show_col_types = FALSE
  ) %>% rename(chr = X1, snp = X2, cM = X3, bp = X4)
flat <- read_tsv(
  glue("{ref.dir}/GRCh38.p14/recomb_maps_1Mb1cM/chr18.map"), 
  col_names = FALSE, show_col_types = FALSE
  ) %>% rename(chr = X1, snp = X2, cM = X3, bp = X4)
maps <- rbind(
  hapmap %>% mutate(map = "hapmap"),
  decode %>% mutate(map = "decode"),
  flat %>% mutate(map = "1Mb1cM")
  ) %>%
  mutate(
    region = case_when(
      bp < cen.start.bp ~ "pre",
      bp > cen.end.bp ~ "post",
      TRUE ~ "cen"
    )
  )

## plot line plot of maps ----
ggplot(maps) +
  geom_line(aes(x = bp, y = cM, color = map), linewidth = 1) +
  scale_color_manual(
    values = c(
      "1Mb1cM" = "darkgrey",
      "decode"       = "darkred",
      "hapmap"  = "darkblue"
    )
  ) +
  ylim(0, 120) +
  xlim(0, 8e7) +
  geom_vline(xintercept = cen.start.bp, linetype = "dashed") +
  geom_vline(xintercept = cen.end.bp, linetype = "dashed") +
  theme_bw(base_size = 28)

## plot dot plot of maps ----
ggplot() +
  # dots for empirical maps
  geom_point(
    data = maps %>% filter(map %in% c("hapmap", "decode")),
    aes(x = bp, y = cM, color = map), size = 1
  ) +
  # line for flat map only
  geom_line(
    data = maps %>% filter(map == "1Mb1cM"),
    aes(x = bp, y = cM, color = map),
    linewidth = 1
  ) +
  scale_color_manual(
    values = c(
      "1Mb1cM" = "darkgrey",
      "decode" = "darkred",
      "hapmap" = "darkblue"
    )
  ) +
  geom_vline(xintercept = cen.start.bp, linetype = "dashed") +
  geom_vline(xintercept = cen.end.bp,   linetype = "dashed") +
  xlim(0, 8e7) +
  ylim(0, 120) +
  theme_bw(base_size = 28)

# slopes per region on maps ----
bin.size <- 1e6
maps.thin <- maps %>%
  mutate(bp.bin = floor(bp / bin.size) * bin.size) %>%
  group_by(map, bp.bin) %>%
  summarise(bp = mean(bp), cM = mean(cM), .groups = "drop") %>%
  mutate(
    region = case_when(
      bp < cen.start.bp ~ "pre",
      bp > cen.end.bp ~ "post",
      TRUE ~ "cen"
    )
  )

maps.slopes <- maps.thin %>%
  arrange(map, bp) %>%
  group_by(map, region) %>%
  mutate(
    d.bp = bp - lag(bp),
    d.cM = cM - lag(cM),
    slope.cM.per.Mb = (d.cM / d.bp) * 1e6
  ) %>%
  ungroup() %>%
  filter(is.finite(slope.cM.per.Mb), d.bp > 0)

aov.pre  <- aov(slope.cM.per.Mb ~ map, data = maps.slopes %>% filter(region == "pre"))
aov.post <- aov(slope.cM.per.Mb ~ map, data = maps.slopes %>% filter(region == "post"))

TukeyHSD(aov.pre)
TukeyHSD(aov.post)

cen.flanks <- maps %>%
  filter(map %in% c("hapmap", "decode")) %>%
  group_by(map) %>%
  summarise(
    bp.left  = max(bp[bp < cen.start.bp]),
    cM.left  = cM[bp == bp.left],
    bp.right = min(bp[bp > cen.end.bp]),
    cM.right = cM[bp == bp.right],
    .groups = "drop"
  )
cen.slopes <- cen.flanks %>%
  mutate(
    delta.cM = cM.right - cM.left,
    delta.bp = bp.right - bp.left,
    slope.cM.per.bp = delta.cM / delta.bp,
    slope.cM.per.Mb = slope.cM.per.bp * 1e6
  )
cen.slopes

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
# 
# # analyze sensitivity performance ----
# sens.focal.pop <- "ADMIX"
# sens.focal.imap <- "hapmap"
# df.sens <- all.stats.per.rep %>%
#   filter(
#     # pop == sens.focal.pop,
#     imap == sens.focal.imap,
#     filt == "ibdends",
#     software == "rapid"
#   ) %>%
#   pivot_longer(
#     cols = matches("^sensitivity(\\.|$)"),
#     names_to = "threshold",
#     values_to = "sensitivity"
#   ) %>%
#   mutate(
#     threshold = case_when(
#       threshold == "sensitivity" ~ "overall",
#       TRUE ~ sub("^sensitivity\\.", "", threshold)
#     ),
#     threshold = factor(threshold, levels = thr.levels)
#   )
# 
# ggplot(df.sens, aes(x = threshold, y = sensitivity)) +
#   geom_boxplot(
#     aes(group = interaction(threshold, software.code), fill = software.code)
#   ) +
#   facet_grid( ~ pop) +
#   scale_x_discrete(labels = thr.labels) +
#   scale_fill_manual(
#     values = c(
#       "h1" = "lightblue1",
#       "h2"  = "dodgerblue",
#       "h3" = "blue",
#       "r1" = "mistyrose",
#       "r2"  = "red",
#       "r3" = "darkred"
#     )
#   ) +
#   labs(
#     x = "IBD length bins (cM)",
#     y = "Sensitivity",
#     # title = glue("Sensitivty of {sens.focal.pop} IBD with a {sens.focal.imap} map"),
#     title = glue("Sensitivty of RaPID + ibd-ends across populations"),
#     fill = NULL
#   ) +
#   theme_bw(base_size = 32) + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
# 
# # analyze ppv performance ----
# ppv.focal.pop <- "ADMIX"
# ppv.focal.imap <- "hapmap"
# df.ppv <- all.stats.per.rep %>%
#   filter(
#     # pop == ppv.focal.pop,
#     imap == ppv.focal.imap,
#     filt == "hot",
#     software == "hapibd"
#   ) %>%
#   pivot_longer(
#     cols = matches("^ppv(\\.|$)"),
#     names_to = "threshold",
#     values_to = "ppv"
#   ) %>%
#   mutate(
#     threshold = case_when(
#       threshold == "ppv" ~ "overall",
#       TRUE ~ sub("^ppv\\.", "", threshold)
#     ),
#     threshold = factor(threshold, levels = thr.levels)
#   )
# 
# ggplot(df.ppv, aes(x = threshold, y = ppv)) +
#   geom_boxplot(
#     aes(group = interaction(threshold, software.code), fill = software.code)
#   ) +
#   facet_grid(~ pop) +
#   scale_x_discrete(labels = thr.labels) +
#   scale_fill_manual(
#     values = c(
#       "h1" = "lightblue1",
#       "h2"  = "dodgerblue",
#       "h3" = "blue",
#       "r1" = "mistyrose",
#       "r2"  = "red",
#       "r3" = "darkred"
#     )
#   ) +
#   labs(
#     x = "IBD length bins (cM)",
#     y = "ppv",
#     # title = glue("ppv of {ppv.focal.pop} IBD with a {ppv.focal.imap} map"),
#     title = glue("ppv of hap-IBD + ibd-ends across populations"),
#     fill = NULL
#   ) +
#   theme_bw(base_size = 32) + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
# 
# # analyze length distribution using ks testing ----
# ## define table of files specific to one pop ----
# len.focal.pop <- "ADMIX"
# imap.all.dir.vec <- names(dir.paths.map)
# file.tbl <- tibble(dir = imap.all.dir.vec) %>%
#   mutate(
#     file = map(dir, ~ list.files(
#       .x, full.names = TRUE,
#       pattern = "_\\d+(\\.(hot|ibdends|ibdendshot))?\\.bed\\.gz$"
#     ))
#   ) %>%
#   unnest(file) %>%
#   filter(!is.na(file))
# file.tbl <- file.tbl %>%
#   filter(
#     !grepl("\\.ibdendshot\\.bed\\.gz$", file),
#     str_detect(file, len.focal.pop))
# 
# ## parallel read + annotate per file ----
# plan(multisession, workers = cores)
# imap.all.calls <- NA
# imap.all.calls <- with_progress({
#   p <- progressor(along = seq_len(nrow(file.tbl)))
#   
#   future_map_dfr(
#     seq_len(nrow(file.tbl)),
#     function(i) {
#       p()
#       
#       dir <- file.tbl$dir[i]
#       f   <- file.tbl$file[i]
#       fname <- basename(f)
#       parts <- strsplit(fname, "_")[[1]]
#       
#       if (parts[1] == "true") {
#         imap <- "true"
#         software <- "true"
#       } else {
#         imap <- parts[3]
#         software <- parts[1]
#       }
#       
#       pop <- parts[length(parts) - 1]
#       rep <- as.integer(sub("(\\d+).*", "\\1", parts[length(parts)]))
#       filt <- sub(
#         ".*_\\d+(?:\\.(hot|ibdends))?\\.bed\\.gz$",
#         "\\1",
#         fname
#       )
#       if (filt == "") filt <- "none"
#       tmp <- readr::read_tsv(
#         f,
#         col_names = TRUE,
#         show_col_types = FALSE,
#         progress = FALSE
#       )
#       
#       tmp$pop <- pop
#       tmp$imap <- imap
#       tmp$filt <- filt
#       tmp$rep <- rep
#       tmp$software <- software
#       tmp$software.code <- unname(dir.paths.map[dir])
#       sc <- tmp$software.code[1]
#       tmp$param.code <- if (sc == "true") "true" else paste0(substr(sc, 1, 1), stringr::str_extract(sc, "\\d+$"))
#       ids <- unrels.tbl$ids[unrels.tbl$rep == rep & unrels.tbl$pop == pop][[1]]
#       tmp$is.unrel <- if (is.null(ids)) FALSE else (tmp$sam1 %in% ids & tmp$sam2 %in% ids)
#       
#       tmp
#     },
#     .options = furrr_options(seed = TRUE)
#   )
# })
# plan(sequential)
# 
# ## ks testing procedure ----
# ks.res <- imap.all.calls %>%   
#   filter(software.code != "true", is.unrel == TRUE, bp.str > cen.end.bp) %>%
#   group_by(pop, rep, software.code, param.code, imap, filt) %>% 
#   summarise(inf = list(cM.len), .groups = "drop") %>%
#   left_join(
#     imap.all.calls %>%
#       filter(software.code == "true") %>%
#       group_by(pop, rep) %>%
#       summarise(tru = list(cM.len), .groups = "drop"),
#     by = c("pop", "rep")
#   ) %>%
#   mutate(
#     p.value = map2_dbl(tru, inf, run.ks),
#     n.inf   = map_int(inf, length),
#     n.true  = map_int(tru, length)
#   ) %>%
#   select(-tru, -inf)
# 
# ks.res <- ks.res %>%
#   group_by(pop) %>%
#   filter(filt != "ibdendshot") %>%
#   mutate(
#     p.adj = p.adjust(p.value, method = "BH"),
#     log.p.adj = log10(p.adj),
#     is.sig = p.adj <= 0.05,
#     filt = factor(filt, levels = c("none", "hot", "ibdends")),
#     imap = factor(imap, levels = c("hapmap", "decode", "1Mb1cM")),
#     software.code = factor(software.code),
#     param.code = factor(param.code),
#   ) %>%
#   ungroup()
# 
# ks.bar <- ks.res %>%
#   mutate(outcome = if_else(p.adj < 0.05, "reject", "don't reject")) %>%
#   count(pop, imap, filt, software.code, param.code, outcome, name = "n") %>%
#   group_by(pop, imap, filt, software.code) %>%
#   mutate(frac = n / sum(n)) %>%
#   ungroup() %>%
#   mutate(
#     sw.family = case_when(
#       str_starts(software.code, "h") ~ "hapibd",
#       str_starts(software.code, "r") ~ "rapid",
#       TRUE ~ NA_character_
#     ),
#     fill.key = interaction(sw.family, outcome, sep = "_")
#   )
# 
# ## plot fraction of null hypo. rejections ----
# ks.bar.ex <- ks.bar %>% filter(filt == "ibdends", imap == "hapmap")
# ggplot(
#   ks.bar,
#   aes(x = param.code, y = frac, fill = fill.key)
# ) +
#   geom_col(position = "stack", width = 0.8, color = "black", linewidth = 0.5) +
#   facet_grid(imap ~ filt) +
#   scale_y_continuous(limits = c(0, 1)) +
#   scale_fill_manual(
#     values = c(
#       "hapibd_don't reject" = "lightblue1",
#       "hapibd_reject"       = "darkblue",
#       "rapid_don't reject"  = "mistyrose",
#       "rapid_reject"        = "darkred"
#     ),
#     labels = c(
#       "hapibd_don't reject" = "hap-IBD: don't reject",
#       "hapibd_reject"       = "hap-IBD: reject",
#       "rapid_don't reject"  = "RaPID: don't reject",
#       "rapid_reject"        = "RaPID: reject"
#     )
#   ) +
#   labs(
#     title = glue(
#       "K-S Test True vs Inferred {len.focal.pop} IBD length distributions (unrelated pairs)"
#     ),
#     x = "Software / parameterization",
#     y = "Fraction of replicates",
#     fill = NULL
#   ) +
#   theme_bw(base_size = 32)
# 
# # plot ibd distributions of rep 1 ----
# ## build table of files for rep 1 ----
# rep1.file.tbl <- tibble(dir = imap.all.dir.vec) %>%
#   mutate(
#     file = map(dir, ~ list.files(
#       .x, full.names = TRUE,
#       pattern = "_1(\\.(hot|ibdends|ibdendshot))?\\.bed\\.gz$"
#     ))
#   ) %>%
#   unnest(file) %>%
#   filter(!is.na(file)) %>%
#   filter(!grepl("\\.ibdendshot\\.bed\\.gz$", file))
# 
# ## parallel read of rep 1 ibd segements ----
# plan(multisession, workers = cores)
# rep1.calls <- NA
# rep1.calls <- with_progress({
#   p <- progressor(along = seq_len(nrow(rep1.file.tbl)))
#   
#   future_map_dfr(
#     seq_len(nrow(rep1.file.tbl)),
#     function(i) {
#       p()
#       
#       dir <- rep1.file.tbl$dir[i]
#       f   <- rep1.file.tbl$file[i]
#       
#       fname <- basename(f)
#       parts <- strsplit(fname, "_")[[1]]
#       
#       if (parts[1] == "true") {
#         imap <- "true"
#         software <- "true"
#       } else {
#         imap <- parts[3]
#         software <- parts[1]
#       }
#       
#       pop <- parts[length(parts) - 1]
#       rep <- 1L  # because we filtered to _1 files
#       
#       filt <- sub(".*_1(?:\\.(hot|ibdends))?\\.bed\\.gz$", "\\1", fname)
#       if (filt == "") filt <- "none"
#       
#       tmp <- readr::read_tsv(
#         f,
#         col_names = TRUE,
#         show_col_types = FALSE,
#         progress = FALSE
#       )
#       
#       tmp$pop <- pop
#       tmp$imap <- imap
#       tmp$filt <- filt
#       tmp$rep <- rep
#       tmp$software <- software
#       tmp$software.code <- unname(dir.paths.map[dir])
#       
#       sc <- tmp$software.code[1]
#       tmp$param.code <- if (sc == "true") "true" else paste0(
#         substr(sc, 1, 1),
#         stringr::str_extract(sc, "\\d+$")
#       )
#       
#       ids <- unrels.tbl$ids[unrels.tbl$rep == rep & unrels.tbl$pop == pop][[1]]
#       tmp$is.unrel <- if (is.null(ids)) FALSE else (tmp$sam1 %in% ids & tmp$sam2 %in% ids)
#       
#       tmp
#     },
#     .options = furrr_options(seed = TRUE)
#   )
# })
# plan(sequential)
# 
# ## plot ibd length distribution ----
# sub.call.focal.pop <- "ADMIX"
# sub.call.y.max <- 250
# sub.calls <- rep1.calls %>%
#   filter(
#     param.code == "true" | 
#       (param.code == "r1" & filt == "ibdends" & imap == "hapmap") | 
#       (param.code == "r1" & filt == "ibdends" & imap == "1Mb1cM"),
#     pop == sub.call.focal.pop,
#     is.unrel == TRUE
#   ) %>%
#   # hard coded facet wrap order
#   mutate(
#     facet = interaction(param.code, imap, filt, sep = ", "),
#     # facet = fct_relevel(
#     #   facet,
#     #   "h1, 1Mb1cM, none",
#     #   "h1, 1Mb1cM, hot",
#     #   "h1, 1Mb1cM, ibdends",
#     #   "h1, hapmap, none",
#     #   "h1, hapmap, hot",
#     #   "h1, hapmap, ibdends",
#     #   "true, true, none",
#     # )
#     facet = fct_relevel(
#       facet,
#       "h1, 1Mb1cM, ibdends",
#       "h1, hapmap, ibdends",
#       "true, true, none",
#       "r1, 1Mb1cM, ibdends",
#       "r1, hapmap, ibdends"
#     )
#   )
# 
# ggplot(sub.calls) +
#   geom_histogram(
#     aes(x = cM.len, fill = param.code), color = "black",  linewidth = 0.5
#   ) + 
#   facet_wrap(~ facet, drop = TRUE) +
#   # coord_cartesian(ylim = c(0, sub.call.y.max)) +
#   scale_fill_manual(
#     values = c(
#       "h1" = "lightblue1",
#       "h2"  = "dodgerblue",
#       "h3" = "blue",
#       "r1" = "mistyrose",
#       "r2"  = "red",
#       "r3" = "darkred"
#       ),
#     labels = NULL
#     ) + 
#   labs(
#     title = glue("Length of {sub.call.focal.pop} IBD segments"),
#     x = "Length (cM)",
#     y = "count",
#     fill = NULL
#   ) +
#   theme_bw(base_size = 32) +
#   theme(legend.position = "none")
# 
# ## plot ibd colored by centromere location ----
# sub.calls.cen <- sub.calls %>%
#   mutate(
#     not.hits.cen = (bp.end < cen.end.bp) | (bp.str > cen.end.bp),
#     hits.cen = if_else(not.hits.cen, "doesn't overlap", "overlaps centromere"),
#     fill.plot = if_else(hits.cen == "overlaps centromere", "cen", param.code),
#     fill.plot = factor(fill.plot, levels = c("true", "h1", "r1", "cen"))
#   )
# 
# ggplot(sub.calls.cen) +
#   geom_histogram(
#     aes(x = cM.len, fill = fill.plot),
#     position = "stack",
#     color = "black",
#     linewidth = 0.5
#   ) +
#   facet_wrap(~ facet, drop = TRUE) +
#   coord_cartesian(ylim = c(0, sub.call.y.max)) +
#   scale_fill_manual(
#     values = c(
#       "h1"   = "lightblue1",
#       "r1"   = "mistyrose",
#       "true" = "lightgrey",
#       "cen"  = "black"
#     )
#   ) +
#   labs(
#     title = glue(
#       "{sub.call.focal.pop} IBD segments overlapping centromere (black bars)"
#     ),
#     x = "Length (cM)",
#     y = "Count",
#     fill = NULL
#   ) +
#   theme_bw(base_size = 32) +
#   theme(legend.position = "none")
# 
# ## plot ibd colored by before centromere ----
# sub.calls.bfcen <- sub.calls %>%
#   mutate(
#     bef.cen = bp.end < cen.start.bp,
#     bef.cen = if_else(bef.cen, "before centromere", "in/after centromere"),
#     fill.plot = if_else(bef.cen == "before centromere", "bfcen", param.code),
#     fill.plot = factor(fill.plot, levels = c("true", "h1", "r1", "bfcen"))
#   )
# 
# ggplot(sub.calls.bfcen) +
#   geom_histogram(
#     aes(x = cM.len, fill = fill.plot),
#     position = "stack",
#     color = "black",
#     linewidth = 0.5
#   ) +
#   facet_wrap(~ facet, drop = TRUE) +
#   coord_cartesian(ylim = c(0, sub.call.y.max)) +
#   scale_fill_manual(
#     values = c(
#       "h1"   = "lightblue1",
#       "r1"   = "mistyrose",
#       "true" = "lightgrey",
#       "bfcen"  = "black"
#     )
#   ) +
#   labs(
#     title = glue(
#       "{sub.call.focal.pop} IBD segments before centromere (black bars)"
#     ),
#     x = "Length (cM)",
#     y = "Count",
#     fill = NULL
#   ) +
#   theme_bw(base_size = 32) +
#   theme(legend.position = "none")
# 
# ## plot ibd colored by after centromere ----
# sub.calls.aftcen <- sub.calls %>%
#   mutate(
#     aft.cen = bp.str > cen.start.bp,
#     aft.cen = if_else(aft.cen, "after centromere", "before/in centromere"),
#     fill.plot = if_else(aft.cen == "after centromere", "aftcen", param.code),
#     fill.plot = factor(fill.plot, levels = c("true", "h1", "r1", "aftcen"))
#   )
# 
# ggplot(sub.calls.aftcen) +
#   geom_histogram(
#     aes(x = cM.len, fill = fill.plot),
#     position = "stack",
#     color = "black",
#     linewidth = 0.5
#   ) +
#   facet_wrap(~ facet, drop = TRUE) +
#   coord_cartesian(ylim = c(0, sub.call.y.max)) +
#   scale_fill_manual(
#     values = c(
#       "h1"   = "lightblue1",
#       "r1"   = "mistyrose",
#       "true" = "lightgrey",
#       "aftcen"  = "black"
#     )
#   ) +
#   labs(
#     title = glue(
#       "{sub.call.focal.pop} IBD segments after centromere (black bars)"
#     ),
#     x = "Length (cM)",
#     y = "Count",
#     fill = NULL
#   ) +
#   theme_bw(base_size = 32) +
#   theme(legend.position = "none")
# 
# # ## enrichment analysis of ibd before, overlapping, and after chromosome ----
# # region.enrich <- imap.all.calls %>%
# #   filter(
# #     is.unrel == TRUE, 
# #     (param.code == "h1" & filt == "hot") | (param.code == "r2" & filt == "ibdends") | param.code == "true"
# #   ) %>%
# #   mutate(
# #     region = factor(
# #       case_when(
# #         bp.end < cen.end.bp ~ "before",
# #         bp.str > cen.start.bp ~ "after",
# #         TRUE               ~ "overlap"
# #       )
# #     )
# #   ) %>%
# #   count(rep, param.code, software, filt, imap, region, name = "n.calls") %>%
# #   group_by(rep, param.code, software, filt, imap) %>%
# #   mutate(
# #     total.calls = sum(n.calls),
# #     prop.calls = n.calls / total.calls,
# #     imap = factor(imap),
# #     param.code = factor(param.code),
# #     filt = factor(filt)
# #   ) %>%
# #   ungroup()
# # 
# # region.lengths <- data.frame(
# #   region = c("before", "overlap", "after"),
# #   length.bp = c(
# #     cen.start.bp,
# #     cen.end.bp - cen.start.bp,
# #     chr18.len - cen.end.bp
# #   )
# # )
# # 
# # region.test.stat <- region.enrich %>%
# #   left_join(region.lengths, by = "region") %>%
# #   mutate(
# #     rate.per.Mb = n.calls / (length.bp / 1e6),
# #     prop.calls = n.calls / total.calls,
# #     region = factor(region, levels = c("before", "overlap", "after"))
# #   )
# # 
# # method.spec <- tibble::tribble(
# #   ~method,  ~param.code, ~filt,
# #   "rapid",  "r2",        "ibdends",
# #   "hapibd", "h1",        "hot"
# # )
# # 
# # true.df <- region.test.stat %>%
# #   filter(param.code == "true") %>%
# #   select(rep, region, true = all_of("prop.calls"))
# # 
# # inf.df <- region.test.stat %>%
# #   inner_join(method.spec, by = c("param.code", "filt")) %>%
# #   filter(imap %in% c("1Mb1cM", "decode", "hapmap")) %>%
# #   select(method, rep, region, imap, inf = all_of("prop.calls"))
# # 
# # delta.df <- inf.df %>%
# #   inner_join(true.df, by = c("rep", "region")) %>%
# #   mutate(delta = inf - true)
# # 
# # tukey.df <- delta.df %>%
# #   group_by(method, region) %>%
# #   group_modify(~{
# #     mod <- aov(delta ~ imap, data = .x)
# #     tk  <- TukeyHSD(mod, "imap")$imap
# #     tibble(
# #       contrast = rownames(tk),
# #       diff = tk[, "diff"],
# #       lwr  = tk[, "lwr"],
# #       upr  = tk[, "upr"],
# #       p.tukey = tk[, "p adj"]
# #     )
# #   }) %>%
# #   ungroup()
# # 
# # tukey.df <- tukey.df %>%
# #   group_by(method, contrast) %>%
# #   mutate(p.adj.regions = p.adjust(p.tukey, method = "BH")) %>%
# #   ungroup()
# # 
# # ggplot(
# #   region.test.stat %>% filter((software == "hapibd" & filt == "hot") | software == "true" & imap == "true"), 
# #   aes(x = region, y = prop.calls, fill = interaction(param.code, imap))) +
# #   geom_boxplot() +
# #   labs(
# #     x = "Location rel. to centromere",
# #     y = "Prop. of IBD segments",
# #     title = glue("Location-specific {len.focal.pop} IBD abundances using hap-IBD (h1, hot)"),
# #     fill = "Maps"
# #   ) +
# #   scale_fill_manual(
# #     values = c(
# #       "h1.1Mb1cM"   = "lightblue1",
# #       "h1.decode"   = "dodgerblue",
# #       "h1.hapmap"   = "darkblue",
# #       "r2.1Mb1cM"   = "mistyrose",
# #       "r2.decode"   = "red",
# #       "r2.hapmap"   = "darkred",
# #       "true.true" = "lightgrey"
# #     ),
# #     labels = c(
# #       "h1.1Mb1cM"   = "flat",
# #       "h1.decode"   = "decode",
# #       "h1.hapmap"   = "hapmap",
# #       "r2.1Mb1cM"   = "flat",
# #       "r2.decode"   = "decode",
# #       "r2.hapmap"   = "hapmap",
# #       "true.true" = "hapmap (truth)"
# #     )
# #   ) +
# #   theme_bw(base_size = 24)
# 
# # analyze ibd score rel true performance ----
# df.reltrue <- all.stats.per.rep %>%
#   filter(
#     imap == "hapmap",
#     filt %in% c("ibdends")
#   ) %>%
#   pivot_longer(
#     cols = matches("^ibd\\.score\\.rel\\.true(\\.|$)"),
#     values_to = "ibd.score.rel.true"
#   )
# 
# ggplot(df.reltrue, aes(x = filt, y = ibd.score.rel.true, color = software.code)) +
#   geom_boxplot(
#     aes(group = interaction(filt, software.code)), fill = "lightgrey"
#   ) +
#   geom_hline(yintercept = 1, linetype = "dashed") +
#   labs(
#     x = "Filtering procedure",
#     y = "Relative IBD score (inferred / true)",
#     color = "Software",
#     title = "Inferred vs True IBD score using ibd-ends"
#   ) +
#   coord_trans(y = "log10") +
#   facet_wrap(~ pop) +
#   theme_bw()
# 
# ## plot true ibd scores ----
# df.true.ibd <- all.stats.per.rep %>%
#   filter(
#     imap == "hapmap",
#     rep == 1,
#     filt == "ibdends",
#     software.code == "h1"
#   ) %>%
#   pivot_longer(
#     cols = matches("^ibd\\.score\\.err(\\.|$)"),
#     values_to = "ibd.score.err"
#   )
# afr.score <- 119695.72
# df.true.ibd <- df.true.ibd %>%
#   mutate(rel.score = true.ibd.score / afr.score)
# ggplot(df.true.ibd, aes(x = pop, y = true.ibd.score)) +
#   geom_col(fill = "grey") +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   labs(
#     x = NULL,
#     y = "True IBD score",
#     title = glue("IBD score"),
#     fill = NULL
#   ) +
#   theme_bw(base_size = 32) 
# ggplot(df.true.ibd, aes(x = pop, y = rel.score)) +
#   geom_col(fill = "grey") +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   labs(
#     x = NULL,
#     y = "True IBD score rel. to AFR",
#     title = glue("IBD score rel. to AFR"),
#     fill = NULL
#   ) +
#   theme_bw(base_size = 32) 
# 
# # analyze ibd score error performance ----
# df.err <- all.stats.per.rep %>%
#   filter(
#     # imap == "hapmap",
#     # filt == "ibdends",
#     imap != "hapmap" | filt != "hot"
#     ) %>%
#   pivot_longer(
#     cols = matches("^ibd\\.score\\.err(\\.|$)"),
#     values_to = "ibd.score.err"
#   )
# df.err.ex <- df.err %>% filter(filt == "ibdends", imap == "hapmap", software.code %in% c("h1", "r1"))
# ggplot(df.err, aes(x = pop, y = ibd.score.err, fill = software.code)) +
#   geom_boxplot() +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   scale_fill_manual(
#     values = c(
#       "h1" = "lightblue1",
#       "h2"  = "dodgerblue",
#       "h3" = "blue",
#       "r1" = "mistyrose",
#       "r2"  = "red",
#       "r3" = "darkred"
#     )
#   ) +
#   labs(
#     x = NULL,
#     y = "IBD score error",
#     title = glue("IBD score error"),
#     fill = NULL
#   ) +
#   facet_grid(filt ~ imap) +
#   # facet_wrap(~ filt + imap) +
#   theme_bw(base_size = 32) 
# 
# # plot rel ibd score error to AFR ----
# 
# df.err <- all.stats.per.rep %>%
#   # keep whatever imap/filt/software you want; don't force hapmap/ibdends if you want facets
#   # filter(imap != "decode", filt != "none", imap != "hapmap" | filt != "hot") %>%
#   filter(imap == "hapmap", filt == "ibdends") %>%
#   group_by(imap, filt, software.code, rep) %>%   # add other ID cols here if needed (e.g., tmap)
#   mutate(
#     # pull the AFR baselines within each (imap,filt,software,rep) group
#     true.afr = true.ibd.score[pop == "AFR"][1],
#     inf.afr  = inf.ibd.score[pop == "AFR"][1],
#     
#     # normalize to AFR (relative-to-AFR scores)
#     true.ibd.score.rel.afr = true.ibd.score / true.afr,
#     inf.ibd.score.rel.afr  = inf.ibd.score  / inf.afr,
#     
#     # error = (true relative) - (inferred relative)
#     ibd.score.rel.afr.err = true.ibd.score.rel.afr - inf.ibd.score.rel.afr
#   ) %>%
#   ungroup() %>%
#   # optional safety: drop groups where AFR is missing or 0 (would create Inf/NA)
#   filter(is.finite(true.ibd.score.rel.afr),
#          is.finite(inf.ibd.score.rel.afr),
#          is.finite(ibd.score.rel.afr.err))
# 
# ggplot(df.err, aes(x = pop, y = ibd.score.rel.afr.err, fill = software.code)) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   scale_fill_manual(
#     values = c(
#       "h1" = "lightblue1",
#       "h2" = "dodgerblue",
#       "h3" = "blue",
#       "r1" = "mistyrose",
#       "r2" = "red",
#       "r3" = "darkred"
#     )
#   ) +
#   labs(
#     x = NULL,
#     y = "Relative IBD score error",
#     title = "IBD score error (relative to AFR)" ,
#     fill = NULL
#   ) +
#   # facet_grid(filt ~ imap) +
#   # facet_wrap(~ filt + imap) +
#   theme_bw(base_size = 32)
# 
# # analyze kindship performance ----
# # ---- set these two dirs correctly ----
# hap.dir <- "~/scratch/ibd_compare_chr18_admix/true_vs_inference/hapibd_sim_tmapHapMapII_GRCh38_imapHapMapII_GRCh38_seed1.0_gap2000_extend1.0_out2.0_mark100_mac2/ibd_kinship_mats"
# rap.dir <- "~/scratch/ibd_compare_chr18_admix/true_vs_inference/rapid_sim_tmapHapMapII_GRCh38_imapHapMapII_GRCh38_d2_w100_r10_s2/ibd_kinship_mats"
# 
# # ---- read mats ----
# true.hap <- readRDS(file.path(hap.dir, "true_AFR_30.hot.rds"))
# inf.hap  <- readRDS(file.path(hap.dir, "inf_AFR_30.hot.rds"))
# true.rap <- readRDS(file.path(rap.dir, "true_AFR_30.hot.rds"))
# inf.rap  <- readRDS(file.path(rap.dir, "inf_AFR_30.hot.rds"))
# 
# # ----------------------------
# # helper: evaluate one method
# # ----------------------------
# eval_two_stage <- function(true_mat, inf_mat, software,
#                            eps = 0,               # threshold for "related" (use >0 or small eps)
#                            ut = NULL) {
#   
#   stopifnot(is.matrix(true_mat) || is.data.frame(true_mat))
#   stopifnot(is.matrix(inf_mat)  || is.data.frame(inf_mat))
#   
#   true_mat <- as.matrix(true_mat)
#   inf_mat  <- as.matrix(inf_mat)
#   
#   if (is.null(ut)) ut <- upper.tri(true_mat, diag = FALSE)
#   
#   true <- as.numeric(true_mat[ut])
#   inf  <- as.numeric(inf_mat[ut])
#   
#   df <- tibble(
#     software = software,
#     true = true,
#     inf  = inf
#   ) %>%
#     filter(is.finite(true), is.finite(inf))
#   
#   # ---- stage (i): detection of related pairs ----
#   df <- df %>%
#     mutate(
#       true_rel = true > eps,
#       inf_rel  = inf  > eps
#     )
#   
#   tp <- as.numeric(sum(df$true_rel & df$inf_rel))
#   fp <- as.numeric(sum(!df$true_rel & df$inf_rel))
#   fn <- as.numeric(sum(df$true_rel & !df$inf_rel))
#   tn <- as.numeric(sum(!df$true_rel & !df$inf_rel))
#   
#   precision <- ifelse(tp + fp == 0, NA_real_, tp / (tp + fp))
#   recall    <- ifelse(tp + fn == 0, NA_real_, tp / (tp + fn))
#   f1        <- ifelse(is.na(precision) | is.na(recall) | (precision + recall) == 0,
#                       NA_real_, 2 * precision * recall / (precision + recall))
#   bal_acc   <- {
#     tpr <- ifelse(tp + fn == 0, NA_real_, tp / (tp + fn))
#     tnr <- ifelse(tn + fp == 0, NA_real_, tn / (tn + fp))
#     ifelse(is.na(tpr) | is.na(tnr), NA_real_, (tpr + tnr) / 2)
#   }
#   mcc <- {
#     tp <- as.numeric(tp); fp <- as.numeric(fp)
#     fn <- as.numeric(fn); tn <- as.numeric(tn)
#     
#     denom <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
#     ifelse(denom == 0, NA_real_, (tp * tn - fp * fn) / denom)
#   }
#   
#   detect_stats <- tibble(
#     software = software,
#     eps = eps,
#     tp = tp, fp = fp, fn = fn, tn = tn,
#     precision = precision,
#     recall = recall,
#     f1 = f1,
#     balanced_accuracy = bal_acc,
#     mcc = mcc,
#     n_pairs = nrow(df),
#     frac_true_related = mean(df$true_rel),
#     frac_inf_related  = mean(df$inf_rel)
#   )
#   
#   # ---- stage (ii): accuracy conditional on relatedness ----
#   # common choice: condition on true-related pairs
#   df_rel <- df %>% filter(true_rel)
#   
#   spearman_r <- ifelse(nrow(df_rel) < 2, NA_real_,
#                        cor(df_rel$true, df_rel$inf, method = "spearman", use = "complete.obs"))
#   
#   # optional: also compute on intersection true_rel & inf_rel
#   df_both_rel <- df %>% filter(true_rel & inf_rel)
#   spearman_r_both <- ifelse(nrow(df_both_rel) < 2, NA_real_,
#                             cor(df_both_rel$true, df_both_rel$inf, method = "spearman", use = "complete.obs"))
#   
#   acc_stats <- tibble(
#     software = software,
#     eps = eps,
#     spearman_r_true_related = spearman_r,
#     n_true_related_pairs = nrow(df_rel),
#     spearman_r_both_related = spearman_r_both,
#     n_both_related_pairs = nrow(df_both_rel)
#   )
#   
#   list(df = df, detect_stats = detect_stats, acc_stats = acc_stats)
# }
# 
# # ----------------------------
# # example: run for hapibd vs rapid
# # ----------------------------
# 
# # assume you already read these:
# # true.hap, inf.hap, true.rap, inf.rap
# 
# eps <- 0  # set to something like 1e-6 if you want a strict "nonzero" threshold
# ut  <- upper.tri(true.hap, diag = FALSE)  # reuse same mask across methods
# 
# res.hap <- eval_two_stage(true.hap, inf.hap, software = "hapibd", eps = eps, ut = ut)
# res.rap <- eval_two_stage(true.rap, inf.rap, software = "rapid",  eps = eps, ut = ut)
# 
# detect.tbl <- bind_rows(res.hap$detect_stats, res.rap$detect_stats)
# acc.tbl    <- bind_rows(res.hap$acc_stats,    res.rap$acc_stats)
# 
# detect.tbl
# acc.tbl
# 
# # ----------------------------
# # plots
# # ----------------------------
# 
# # (i) detection: precision/recall/F1 per software
# detect.long <- detect.tbl %>%
#   select(software, precision, recall, f1, balanced_accuracy, mcc) %>%
#   pivot_longer(-software, names_to = "metric", values_to = "value")
# 
# ggplot(detect.long, aes(x = metric, y = value)) +
#   geom_col() +
#   facet_wrap(~ software) +
#   coord_cartesian(ylim = c(0, 1)) +
#   theme_bw() +
#   labs(title = "Stage (i): detection of related pairs", x = NULL, y = NULL)
# 
# # (ii) accuracy conditional on relatedness: scatter on true-related pairs + Spearman
# df_plot <- bind_rows(
#   res.hap$df %>% mutate(software = "hapibd"),
#   res.rap$df %>% mutate(software = "rapid")
# ) %>%
#   mutate(true_rel = true > eps) %>%
#   filter(true_rel)
# 
# lab_plot <- acc.tbl %>%
#   transmute(
#     software,
#     label = sprintf("Spearman r = %.3f\n(n = %d)", spearman_r_true_related, n_true_related_pairs)
#   )
# 
# ggplot(df_plot, aes(x = true, y = inf)) +
#   geom_point(alpha = 0.35, size = 1) +
#   facet_wrap(~ software, scales = "free") +
#   geom_text(
#     data = lab_plot,
#     aes(x = -Inf, y = Inf, label = label),
#     inherit.aes = FALSE,
#     hjust = -0.05, vjust = 1.2
#   ) +
#   theme_bw() +
#   labs(
#     title = "Stage (ii): accuracy conditional on true relatedness",
#     x = "True kinship",
#     y = "Inferred kinship"
#   )
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # test this ----
# 
# 
# 
# # inputs (your dirs) ----
# hap.hap.dir <- "/home1/karatas/scratch/ibd_compare_chr18_admix/true_vs_inference/hapibd_sim_tmapHapMapII_GRCh38_imapHapMapII_GRCh38_seed1.0_gap2000_extend1.0_out2.0_mark100_mac2/ibd_kinship_mats"
# hap.fla.dir <- "/home1/karatas/scratch/ibd_compare_chr18_admix/true_vs_inference/hapibd_sim_tmapHapMapII_GRCh38_imapFlat1Mb1cM_seed1.0_gap2000_extend1.0_out2.0_mark100_mac2/ibd_kinship_mats"
# rap.hap.dir <- "/home1/karatas/scratch/ibd_compare_chr18_admix/true_vs_inference/rapid_sim_tmapHapMapII_GRCh38_imapHapMapII_GRCh38_d2_w100_r10_s2/ibd_kinship_mats"
# rap.fla.dir <- "/home1/karatas/scratch/ibd_compare_chr18_admix/true_vs_inference/rapid_sim_tmapHapMapII_GRCh38_imapFlat1Mb1cM_d2_w100_r10_s2/ibd_kinship_mats"
# 
# kin.dirs <- c(hap.hap.dir, hap.fla.dir, rap.hap.dir, rap.fla.dir)
# 
# # settings ----
# cores <- 8
# use.only.nonzero <- FALSE  # TRUE => keep only pairs with true>0 & inf>0
# use.log10 <- FALSE         # TRUE => cor(log10(.)) on positives only
# 
# # helper: parse metadata from parent directory name ----
# parse.kin.dir.meta <- function(kin.dir) {
#   parent <- basename(dirname(kin.dir))
#   
#   software <- case_when(
#     str_detect(parent, "^hapibd_") ~ "hapibd",
#     str_detect(parent, "^rapid_")  ~ "rapid",
#     TRUE                           ~ "unknown"
#   )
#   
#   tmap <- str_match(parent, "tmap([^_]+)")[, 2]
#   imap <- str_match(parent, "imap([^_]+)")[, 2]
#   
#   params <- NA_character_
#   if (!is.na(imap)) {
#     params <- str_replace(parent, paste0("^.*_imap", imap, "_"), "")
#   }
#   
#   tibble(
#     kin.dir = kin.dir,
#     parent.dir = parent,
#     software = software,
#     tmap = tolower(coalesce(tmap, NA_character_)),
#     imap = tolower(coalesce(imap, NA_character_)),
#     software.with.params = paste0(software, "_", params)
#   )
# }
# 
# # helper: build paired true/inf table for one dir ----
# build.kin.file.table <- function(kin.dir) {
#   files <- list.files(kin.dir, pattern = "\\.rds$", full.names = TRUE)
#   if (length(files) == 0) return(tibble())
#   
#   tbl <- tibble(file = files, fname = basename(files)) %>%
#     mutate(
#       type = str_match(fname, "^(true|inf)_")[, 2],
#       pop  = str_match(fname, "^(?:true|inf)_([^_]+)_")[, 2],
#       rep  = as.integer(str_match(fname, "^(?:true|inf)_[^_]+_(\\d+)")[, 2]),
#       filt = str_match(fname, "^(?:true|inf)_[^_]+_\\d+\\.([^.]+)\\.rds$")[, 2]
#     ) %>%
#     filter(!is.na(type), !is.na(pop), !is.na(rep), !is.na(filt)) %>%
#     select(type, pop, rep, filt, file)
#   
#   tbl %>%
#     pivot_wider(names_from = type, values_from = file) %>%
#     filter(!is.na(true), !is.na(inf)) %>%
#     arrange(pop, filt, rep)
# }
# 
# # helper: compute correlations for one pair ----
# compute.kin.cors.one <- function(true.file, inf.file,
#                                  use.only.nonzero = FALSE,
#                                  use.log10 = FALSE) {
#   k.true <- readRDS(true.file)
#   k.inf  <- readRDS(inf.file)
#   
#   if (!is.matrix(k.true) || !is.matrix(k.inf)) {
#     stop("One of the RDS files is not a matrix: ", true.file, " / ", inf.file)
#   }
#   if (!all(dim(k.true) == dim(k.inf))) {
#     stop("Matrix dim mismatch: ", true.file, " vs ", inf.file)
#   }
#   
#   ut <- upper.tri(k.true, diag = FALSE)
#   true.vec <- k.true[ut]
#   inf.vec  <- k.inf[ut]
#   
#   keep <- is.finite(true.vec) & is.finite(inf.vec)
#   if (use.only.nonzero) keep <- keep & (true.vec > 0) & (inf.vec > 0)
#   
#   true.k <- true.vec[keep]
#   inf.k  <- inf.vec[keep]
#   
#   if (length(true.k) < 3L) {
#     return(tibble(n.pairs = length(true.k), cor.raw = NA_real_, cor.proc = NA_real_))
#   }
#   
#   cor.raw <- suppressWarnings(cor(true.k, inf.k, use = "complete.obs"))
#   
#   if (use.log10) {
#     keep.pos <- (true.k > 0) & (inf.k > 0)
#     true.l <- log10(true.k[keep.pos])
#     inf.l  <- log10(inf.k[keep.pos])
#     cor.proc <- if (length(true.l) < 3L) NA_real_ else suppressWarnings(cor(true.l, inf.l))
#   } else {
#     cor.proc <- cor.raw
#   }
#   
#   tibble(n.pairs = length(true.k), cor.raw = cor.raw, cor.proc = cor.proc)
# }
# 
# # --------------------------------------------------------------------------
# # main code body (no wrapper function)
# # --------------------------------------------------------------------------
# 
# handlers(global = TRUE)
# plan(multisession, workers = cores)
# 
# dir.meta <- map_dfr(kin.dirs, parse.kin.dir.meta)
# 
# file.tbl <- map_dfr(kin.dirs, function(d) {
#   ft <- build.kin.file.table(d)
#   if (nrow(ft) == 0) return(tibble())
#   ft %>% mutate(kin.dir = d)
# }) %>%
#   left_join(dir.meta, by = "kin.dir")
# 
# if (nrow(file.tbl) == 0) {
#   plan(sequential)
#   stop("No paired true/inf RDS files found. (Do you have matching inf_* files?)")
# }
# 
# message("paired true/inf files: ", nrow(file.tbl))
# 
# cors.df <- with_progress({
#   p <- progressor(along = seq_len(nrow(file.tbl)))
#   
#   future_map_dfr(
#     seq_len(nrow(file.tbl)),
#     function(i) {
#       p()
#       
#       row <- file.tbl[i, ]
#       
#       cors <- compute.kin.cors.one(
#         true.file = row$true,
#         inf.file  = row$inf,
#         use.only.nonzero = use.only.nonzero,
#         use.log10 = use.log10
#       )
#       
#       bind_cols(
#         row %>% select(
#           software, software.with.params, tmap, imap,
#           pop, rep, filt, parent.dir, kin.dir
#         ),
#         cors
#       )
#     },
#     .options = furrr_options(seed = TRUE)
#   )
# })
# 
# plan(sequential)
# 
# # quick check
# dplyr::glimpse(cors.df)
# 
# # plotting examples ----
# library(ggplot2)
# 
# ggplot(cors.df, aes(x = software.with.params, y = cor.proc)) +
#   geom_boxplot(outlier.size = 0.6) +
#   geom_jitter(width = 0.15, height = 0, alpha = 0.35, size = 0.6) +
#   facet_grid(pop ~ filt, scales = "free_y") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   labs(
#     x = NULL,
#     y = "Correlation (upper-tri, no diag)",
#     title = "Kinship correlation across replicates"
#   )
# 
# ggplot(cors.df, aes(x = rep, y = cor.proc, group = interaction(pop, filt))) +
#   geom_line(alpha = 0.3) +
#   geom_point(size = 0.7, alpha = 0.5) +
#   facet_grid(pop ~ software.with.params) +
#   theme_bw() +
#   labs(
#     x = "Replicate",
#     y = "Correlation",
#     title = "Per-replicate kinship correlation (true vs inferred)"
#   )
