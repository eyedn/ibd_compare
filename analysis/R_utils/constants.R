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

n.samples <- 300

# chr18
cen.start.bp <- 15460900
cen.end.bp <- 20861206
chr18.len <- 80373285
chr18.map.len <- data.frame(
  "hapmap" = c(117.708923, 80257297),
  "flat"   = c(92.618155, 80373285),
  "decode" = c(105.783407, 80257297)
)

# various directories frequency used
ibd.comp.dir <- "/home1/karatas/scratch/ibd_compare_chr18_admix"
ref.dir <- "/home1/karatas/references"
sim.dir <- "/home1/karatas/scratch/stdpopsim_runs_chr18_admix"

# all stats per replicate column names
stats.meta.names <- c(
  "software", "software.with.params", "tmap", "imap", "filt", "pop", "rep",
  "n.true.unrels", "n.inf.unrels", "true.ibd.score", "inf.ibd.score",
  "ibd.score.rel.true", "ibd.score.err"
)

# for plotting where the x axis is composed of discrete ibd length bins ----
thr.levels <- c("overall", "2", "2.5", "3", "3.5", "4", "4.5", "5")
thr.labels <- c(
  "overall" = "overall",
  "2"   = "2–2.5",
  "2.5" = "2.5–3",
  "3"   = "3–3.5",
  "3.5" = "3.5–4",
  "4"   = "4–4.5",
  "4.5" = "4.5–5",
  "5"   = "≥5"
)

thr.labels.descriptive <- c(
  "overall" = "All IBD",
  "2"   = "Short IBD\n(2–2.5 cM)",
  "5"   = "Long IBD\n(≥5 cM)"
)

# map software_param string to software.code
param.map <- c(
  "hapibd_seed1.0_gap1000_extend1.0_out2.0_mark100_mac2" = "h1",
  "hapibd_seed1.0_gap2000_extend1.0_out2.0_mark100_mac2" = "h2",
  "hapibd_seed2.0_gap1000_extend1.0_out2.0_mark100_mac2" = "h3",
  "rapid_d2_w100_r10_s2" = "r1",
  "rapid_d2_w200_r10_s2" = "r2",
  "rapid_d2_w300_r10_s2" = "r3",
  "true" = "true"
)

# map sub-dirs under ibd.comp.dir to software.code
dir.paths.map <- c(
  setNames("hh1", glue("{ibd.comp.dir}/hapibd_sim_tmapHapMapII_GRCh38_imapHapMapII_GRCh38_seed1.0_gap1000_extend1.0_out2.0_mark100_mac2")),
  setNames("hh2", glue("{ibd.comp.dir}/hapibd_sim_tmapHapMapII_GRCh38_imapHapMapII_GRCh38_seed1.0_gap2000_extend1.0_out2.0_mark100_mac2")),
  setNames("hh3", glue("{ibd.comp.dir}/hapibd_sim_tmapHapMapII_GRCh38_imapHapMapII_GRCh38_seed2.0_gap1000_extend1.0_out2.0_mark100_mac2")),

  setNames("hd1", glue("{ibd.comp.dir}/hapibd_sim_tmapHapMapII_GRCh38_imapDeCodeSexAveraged_GRCh38_seed1.0_gap1000_extend1.0_out2.0_mark100_mac2")),
  setNames("hd2", glue("{ibd.comp.dir}/hapibd_sim_tmapHapMapII_GRCh38_imapDeCodeSexAveraged_GRCh38_seed1.0_gap2000_extend1.0_out2.0_mark100_mac2")),
  setNames("hd3", glue("{ibd.comp.dir}/hapibd_sim_tmapHapMapII_GRCh38_imapDeCodeSexAveraged_GRCh38_seed2.0_gap1000_extend1.0_out2.0_mark100_mac2")),
  
  setNames("hf1", glue("{ibd.comp.dir}/hapibd_sim_tmapHapMapII_GRCh38_imapFlat1Mb1cM_seed1.0_gap1000_extend1.0_out2.0_mark100_mac2")),
  setNames("hf2", glue("{ibd.comp.dir}/hapibd_sim_tmapHapMapII_GRCh38_imapFlat1Mb1cM_seed1.0_gap2000_extend1.0_out2.0_mark100_mac2")),
  setNames("hf3", glue("{ibd.comp.dir}/hapibd_sim_tmapHapMapII_GRCh38_imapFlat1Mb1cM_seed2.0_gap1000_extend1.0_out2.0_mark100_mac2")),
  
  setNames("rh1", glue("{ibd.comp.dir}/rapid_sim_tmapHapMapII_GRCh38_imapHapMapII_GRCh38_d2_w100_r10_s2")),
  setNames("rh2", glue("{ibd.comp.dir}/rapid_sim_tmapHapMapII_GRCh38_imapHapMapII_GRCh38_d2_w200_r10_s2")),
  setNames("rh3", glue("{ibd.comp.dir}/rapid_sim_tmapHapMapII_GRCh38_imapHapMapII_GRCh38_d2_w300_r10_s2")),
  
  setNames("rd1", glue("{ibd.comp.dir}/rapid_sim_tmapHapMapII_GRCh38_imapDeCodeSexAveraged_GRCh38_d2_w100_r10_s2")),
  setNames("rd2", glue("{ibd.comp.dir}/rapid_sim_tmapHapMapII_GRCh38_imapDeCodeSexAveraged_GRCh38_d2_w200_r10_s2")),
  setNames("rd3", glue("{ibd.comp.dir}/rapid_sim_tmapHapMapII_GRCh38_imapDeCodeSexAveraged_GRCh38_d2_w300_r10_s2")),
  
  setNames("rf1", glue("{ibd.comp.dir}/rapid_sim_tmapHapMapII_GRCh38_imapFlat1Mb1cM_d2_w100_r10_s2")),
  setNames("rf2", glue("{ibd.comp.dir}/rapid_sim_tmapHapMapII_GRCh38_imapFlat1Mb1cM_d2_w200_r10_s2")),
  setNames("rf3", glue("{ibd.comp.dir}/rapid_sim_tmapHapMapII_GRCh38_imapFlat1Mb1cM_d2_w300_r10_s2")),
  
  setNames("true", glue("{ibd.comp.dir}/true_ibd"))
)

# map sub-dirs under ibd.comp.dir/true_vs_inference to software.code
tvi.dir.paths.map <- c(
  setNames("hh1", glue("{ibd.comp.dir}/true_vs_inference/hapibd_sim_tmapHapMapII_GRCh38_imapHapMapII_GRCh38_seed1.0_gap1000_extend1.0_out2.0_mark100_mac2")),
  setNames("hh2", glue("{ibd.comp.dir}/true_vs_inference/hapibd_sim_tmapHapMapII_GRCh38_imapHapMapII_GRCh38_seed1.0_gap2000_extend1.0_out2.0_mark100_mac2")),
  setNames("hh3", glue("{ibd.comp.dir}/true_vs_inference/hapibd_sim_tmapHapMapII_GRCh38_imapHapMapII_GRCh38_seed2.0_gap1000_extend1.0_out2.0_mark100_mac2")),
  
  setNames("hd1", glue("{ibd.comp.dir}/true_vs_inference/hapibd_sim_tmapHapMapII_GRCh38_imapDeCodeSexAveraged_GRCh38_seed1.0_gap1000_extend1.0_out2.0_mark100_mac2")),
  setNames("hd2", glue("{ibd.comp.dir}/true_vs_inference/hapibd_sim_tmapHapMapII_GRCh38_imapDeCodeSexAveraged_GRCh38_seed1.0_gap2000_extend1.0_out2.0_mark100_mac2")),
  setNames("hd3", glue("{ibd.comp.dir}/true_vs_inference/hapibd_sim_tmapHapMapII_GRCh38_imapDeCodeSexAveraged_GRCh38_seed2.0_gap1000_extend1.0_out2.0_mark100_mac2")),
  
  setNames("hf1", glue("{ibd.comp.dir}/true_vs_inference/hapibd_sim_tmapHapMapII_GRCh38_imapFlat1Mb1cM_seed1.0_gap1000_extend1.0_out2.0_mark100_mac2")),
  setNames("hf2", glue("{ibd.comp.dir}/true_vs_inference/hapibd_sim_tmapHapMapII_GRCh38_imapFlat1Mb1cM_seed1.0_gap2000_extend1.0_out2.0_mark100_mac2")),
  setNames("hf3", glue("{ibd.comp.dir}/true_vs_inference/hapibd_sim_tmapHapMapII_GRCh38_imapFlat1Mb1cM_seed2.0_gap1000_extend1.0_out2.0_mark100_mac2")),
  
  setNames("rh1", glue("{ibd.comp.dir}/true_vs_inference/rapid_sim_tmapHapMapII_GRCh38_imapHapMapII_GRCh38_d2_w100_r10_s2")),
  setNames("rh2", glue("{ibd.comp.dir}/true_vs_inference/rapid_sim_tmapHapMapII_GRCh38_imapHapMapII_GRCh38_d2_w200_r10_s2")),
  setNames("rh3", glue("{ibd.comp.dir}/true_vs_inference/rapid_sim_tmapHapMapII_GRCh38_imapHapMapII_GRCh38_d2_w300_r10_s2")),
  
  setNames("rd1", glue("{ibd.comp.dir}/true_vs_inference/rapid_sim_tmapHapMapII_GRCh38_imapDeCodeSexAveraged_GRCh38_d2_w100_r10_s2")),
  setNames("rd2", glue("{ibd.comp.dir}/true_vs_inference/rapid_sim_tmapHapMapII_GRCh38_imapDeCodeSexAveraged_GRCh38_d2_w200_r10_s2")),
  setNames("rd3", glue("{ibd.comp.dir}/true_vs_inference/rapid_sim_tmapHapMapII_GRCh38_imapDeCodeSexAveraged_GRCh38_d2_w300_r10_s2")),
  
  setNames("rf1", glue("{ibd.comp.dir}/true_vs_inference/rapid_sim_tmapHapMapII_GRCh38_imapFlat1Mb1cM_d2_w100_r10_s2")),
  setNames("rf2", glue("{ibd.comp.dir}/true_vs_inference/rapid_sim_tmapHapMapII_GRCh38_imapFlat1Mb1cM_d2_w200_r10_s2")),
  setNames("rf3", glue("{ibd.comp.dir}/true_vs_inference/rapid_sim_tmapHapMapII_GRCh38_imapFlat1Mb1cM_d2_w300_r10_s2"))
)

# map sub-dirs where imap == hapmap under ibd.comp.dir to software.code ----
imap.hapmap.dir.paths.map <- c(
  setNames("hh1", glue("{ibd.comp.dir}/hapibd_sim_tmapHapMapII_GRCh38_imapHapMapII_GRCh38_seed1.0_gap1000_extend1.0_out2.0_mark100_mac2")),
  setNames("hh2", glue("{ibd.comp.dir}/hapibd_sim_tmapHapMapII_GRCh38_imapHapMapII_GRCh38_seed1.0_gap2000_extend1.0_out2.0_mark100_mac2")),
  setNames("hh3", glue("{ibd.comp.dir}/hapibd_sim_tmapHapMapII_GRCh38_imapHapMapII_GRCh38_seed2.0_gap1000_extend1.0_out2.0_mark100_mac2")),
  
  setNames("rh1", glue("{ibd.comp.dir}/rapid_sim_tmapHapMapII_GRCh38_imapHapMapII_GRCh38_d2_w100_r10_s2")),
  setNames("rh2", glue("{ibd.comp.dir}/rapid_sim_tmapHapMapII_GRCh38_imapHapMapII_GRCh38_d2_w200_r10_s2")),
  setNames("rh3", glue("{ibd.comp.dir}/rapid_sim_tmapHapMapII_GRCh38_imapHapMapII_GRCh38_d2_w300_r10_s2")),
  
  setNames("true", glue("{ibd.comp.dir}/true_ibd"))
)