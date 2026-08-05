# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# constants.R
# ______________________________________________________________________________

# pattern: Functional Core

library(glue)

n.samples <- 500
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

ibd.comp.dir <- "/home1/karatas/scratch/ibd_compare_chr18_ooa/filter_format"
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

software.palette <- c("h3" = "blue", "r1" = "mistyrose")
pop.palette <- c("ADX" = "dodgerblue", "AFR" = "mistyrose", "EUR" = "darkred")
filt.levels <- c("none", "hot", "ibdends", "hotibdends", "ibdendshot")
filt.labels <- c("none" = "N", "hot" = "H", "ibdends" = "I",
                 "hotibdends" = "HI", "ibdendshot" = "IH")
filt.labeller <- ggplot2::labeller(filt = filt.labels)

param.map <- c(
  "hapibd_minseed2.0_gap1000_extend1.0_out2.0_mark100_mac2" = "h3",
  "rapid_d2_w100_r10_s2" = "r1", "true" = "true"
)

dir.paths.map <- c(
  setNames("hh3", glue("{ibd.comp.dir}/hapibd_sim_tmapHapMapII_GRCh38_imapHapMapII_GRCh38_minseed2.0_gap1000_extend1.0_out2.0_mark100_mac2")),
  setNames("hd3", glue("{ibd.comp.dir}/hapibd_sim_tmapHapMapII_GRCh38_imapDeCodeSexAveraged_GRCh38_minseed2.0_gap1000_extend1.0_out2.0_mark100_mac2")),
  setNames("hf3", glue("{ibd.comp.dir}/hapibd_sim_tmapHapMapII_GRCh38_imapFlat1Mb1cM_minseed2.0_gap1000_extend1.0_out2.0_mark100_mac2")),
  setNames("rh1", glue("{ibd.comp.dir}/rapid_sim_tmapHapMapII_GRCh38_imapHapMapII_GRCh38_d2_w100_r10_s2")),
  setNames("rd1", glue("{ibd.comp.dir}/rapid_sim_tmapHapMapII_GRCh38_imapDeCodeSexAveraged_GRCh38_d2_w100_r10_s2")),
  setNames("rf1", glue("{ibd.comp.dir}/rapid_sim_tmapHapMapII_GRCh38_imapFlat1Mb1cM_d2_w100_r10_s2")),
  setNames("true", glue("{ibd.comp.dir}/true_ibd"))
)

tvi.dir.paths.map <- c(
  setNames("hh3", glue("{ibd.comp.dir}/true_vs_inference/hapibd_sim_tmapHapMapII_GRCh38_imapHapMapII_GRCh38_minseed2.0_gap1000_extend1.0_out2.0_mark100_mac2")),
  setNames("hf3", glue("{ibd.comp.dir}/true_vs_inference/hapibd_sim_tmapHapMapII_GRCh38_imapFlat1Mb1cM_minseed2.0_gap1000_extend1.0_out2.0_mark100_mac2")),
  setNames("rh1", glue("{ibd.comp.dir}/true_vs_inference/rapid_sim_tmapHapMapII_GRCh38_imapHapMapII_GRCh38_d2_w100_r10_s2")),
  setNames("rf1", glue("{ibd.comp.dir}/true_vs_inference/rapid_sim_tmapHapMapII_GRCh38_imapFlat1Mb1cM_d2_w100_r10_s2"))
)

imap.hapmap.dir.paths.map <- c(
  setNames("hh3", glue("{ibd.comp.dir}/hapibd_sim_tmapHapMapII_GRCh38_imapHapMapII_GRCh38_minseed2.0_gap1000_extend1.0_out2.0_mark100_mac2")),
  setNames("rh1", glue("{ibd.comp.dir}/rapid_sim_tmapHapMapII_GRCh38_imapHapMapII_GRCh38_d2_w100_r10_s2")),
  setNames("true", glue("{ibd.comp.dir}/true_ibd"))
)
