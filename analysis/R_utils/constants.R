# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# constants.R
# ______________________________________________________________________________


library(glue)

SIM_DIR <- Sys.getenv(
  "IBD_COMPARE_SIM_DIR", "/home1/karatas/scratch/ibd_compare_chr18_ooa"
)
POP_INFO_DIR <- file.path(SIM_DIR, "pop_info")
CHROMOSOME <- 18L
ANALYSIS_SAMPLE_SIZE <- 50L
FOCAL_REPLICATE <- 1L
ANALYSIS_POPULATIONS <- c("ADX", "AFR", "EUR")
WORKERS <- Sys.getenv("SLURM_CPUS_PER_TASK")
WORKERS <- if (grepl("^[1-9][0-9]*$", WORKERS)) {
  suppressWarnings(as.integer(WORKERS))
} else 1L
if (is.na(WORKERS) || WORKERS < 1L) WORKERS <- 1L
sim.dir <- SIM_DIR
pop.info.dir <- POP_INFO_DIR
n.samples <- ANALYSIS_SAMPLE_SIZE
ibd.score.mse.epsilon <- 1e-8
kinship.mse.epsilon <- 1e-8
kinship.corr.clamp.epsilon <- 1e-8

cen.start.bp <- 15460900
cen.end.bp <- 20861206
chr18.len <- 80373285
chr18.map.len <- data.frame(
  "hapmap" = c(117.708923, 80257297),
  "flat" = c(92.618155, 80373285),
  "decode" = c(105.783407, 80257297)
)

ibd.comp.dir <- file.path(SIM_DIR, "filter_format")
ref.dir <- "/home1/karatas/references"

stats.meta.names <- c(
  "software", "software.with.params", "tmap", "imap", "filt", "pop", "rep",
  "n.true.unrels", "n.inf.unrels", "true.ibd.score", "inf.ibd.score",
  "ibd.score.rel.true", "ibd.score.err"
)

thr.levels <- c(
  "overall", "2", "2.5", "3", "3.5", "4", "4.5", "5", "5.5",
  "6", "6.5", "7", "7.5", "8", "8.5", "9", "9.5", "10"
)
thr.labels <- c(
  "overall" = "overall", "2" = "2-2.5", "2.5" = "2.5-3", "3" = "3-3.5",
  "3.5" = "3.5-4", "4" = "4-4.5", "4.5" = "4.5-5", "5" = "5-5.5",
  "5.5" = "5.5-6", "6" = "6-6.5", "6.5" = "6.5-7", "7" = "7-7.5",
  "7.5" = "7.5-8", "8" = "8-8.5", "8.5" = "8.5-9", "9" = "9-9.5",
  "9.5" = "9.5-10", "10" = ">10"
)
thr.labels.descriptive <- c(
  "overall" = "All IBD", "2" = "Short IBD\n(2-2.5 cM)",
  "10" = "Long IBD\n(>10 cM)"
)

software.palette <- c(
  h1 = "lightblue1", h2 = "dodgerblue", h3 = "blue",
  r1 = "mistyrose", r2 = "red", r3 = "darkred"
)
pop.palette <- c("ADX" = "dodgerblue", "AFR" = "mistyrose", "EUR" = "darkred")
filt.levels <- c("none", "hot", "ibdends", "hotibdends", "ibdendshot")
filt.labels <- c("none" = "N", "hot" = "H", "ibdends" = "I",
                 "hotibdends" = "HI", "ibdendshot" = "IH")
filt.labeller <- ggplot2::labeller(filt = filt.labels)

# One profile table supplies all parameter and directory mappings.
ANALYSIS_PROFILES <- data.frame(
  code = c("h1", "h2", "h3", "r1", "r2", "r3"),
  profile = c("seed1_gap1000", "seed1_gap2000", "seed2_gap1000",
              "w100", "w200", "w300"),
  software = c(rep("hapibd", 3), rep("rapid", 3)),
  parameters = c(
    "minseed1.0_gap1000_extend1.0_out2.0_mark100_mac2",
    "minseed1.0_gap2000_extend1.0_out2.0_mark100_mac2",
    "minseed2.0_gap1000_extend1.0_out2.0_mark100_mac2",
    "d2_w100_r10_s2", "d2_w200_r10_s2", "d2_w300_r10_s2"
  )
)
ANALYSIS_PROFILES$software.with.params <- with(
  ANALYSIS_PROFILES, paste(software, parameters, sep = "_")
)
INFERENCE_MAPS <- c(
  h = "HapMapII_GRCh38", d = "DeCodeSexAveraged_GRCh38", f = "Flat1Mb1cM"
)
param.map <- c(setNames(ANALYSIS_PROFILES$code,
                       ANALYSIS_PROFILES$software.with.params), true = "true")
PROFILE_DIRS <- unlist(lapply(seq_len(nrow(ANALYSIS_PROFILES)), function(i) {
  profile <- ANALYSIS_PROFILES[i, ]
  dirs <- paste0(profile$software, "_sim_tmapHapMapII_GRCh38_imap",
                 INFERENCE_MAPS, "_", profile$parameters)
  codes <- paste0(substr(profile$code, 1, 1), names(INFERENCE_MAPS),
                  substr(profile$code, 2, 2))
  setNames(codes, dirs)
}), use.names = TRUE)
dir.paths.map <- c(
  setNames(unname(PROFILE_DIRS), file.path(ibd.comp.dir, names(PROFILE_DIRS))),
  setNames("true", file.path(ibd.comp.dir, "true_ibd"))
)
tvi.dir.paths.map <- setNames(
  unname(PROFILE_DIRS),
  file.path(ibd.comp.dir, "true_vs_inference", names(PROFILE_DIRS))
)
imap.hapmap.dir.paths.map <- dir.paths.map[
  dir.paths.map == "true" | substr(dir.paths.map, 2, 2) == "h"
]
PROFILE_PAIRINGS <- expand.grid(
  hapibd = ANALYSIS_PROFILES$code[1:3],
  rapid = ANALYSIS_PROFILES$code[4:6], stringsAsFactors = FALSE
)
rownames(PROFILE_PAIRINGS) <- with(
  PROFILE_PAIRINGS, paste(hapibd, rapid, sep = "_")
)
