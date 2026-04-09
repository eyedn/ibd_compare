# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology 
# Mooney Lab
# ___
# conf_mat.R
# ______________________________________________________________________________


library(tidyverse)
library(glue)
library(lme4)
library(emmeans)

# source functions ----
source("R_utils/constants.R")
source("R_utils/log.msg.R")
source("R_utils/read.all.stats.per.rep.R")
source("R_utils/pivot.confusion.metric.long.R")
source("R_utils/build.confusion.model.formula.R")
source("R_utils/build.confusion.test.specs.R")
source("R_utils/fit.confusion.glmm.R")
source("R_utils/run.confusion.lrt.R")
source("R_utils/tidy.confusion.emmeans.R")
source("R_utils/tidy.confusion.odds.ratios.R")
source("R_utils/run.overall.confusion.analysis.R")
source("R_utils/run.bin.effects.confusion.analysis.R")
source("R_utils/plot.confusion.boxplots.R")

# define constants ----
adj.method <- "bonferroni"
per.rep.path <- glue("{ibd.comp.dir}/true_vs_inference/all_stats.per_rep.tsv")
software.palette <- c(
  "h1" = "lightblue1",
  "h2" = "dodgerblue",
  "h3" = "blue",
  "r1" = "mistyrose",
  "r2" = "red",
  "r3" = "darkred"
)
pop.palette <- c(
  "ADMIX" = "dodgerblue",
  "AFR" = "mistyrose",
  "ASIA" = "red",
  "EUR" = "darkred"
)

# read in replicate-level statistics ----
all.stats.per.rep <- read.all.stats.per.rep(
  per.rep.path,
  stats.meta.names, param.map
)

# build long-format metric tables ----
sens.long <- pivot.confusion.metric.long(
  all.stats.per.rep = all.stats.per.rep,
  metric = "sensitivity",
  thr.levels = thr.levels
  )
precision.long <- pivot.confusion.metric.long(
  all.stats.per.rep = all.stats.per.rep,
  metric = "precision",
  thr.levels = thr.levels
  )

sens.long.overall <- sens.long %>% 
  filter(threshold == "overall", filt != "ibdendshot")
precision.long.overall <- precision.long %>% 
  filter(threshold == "overall", filt != "ibdendshot")

sens.long.bins <- sens.long %>% 
  filter(threshold != "overall", filt != "ibdendshot")
precision.long.bins <- precision.long %>% 
  filter(threshold != "overall", filt != "ibdendshot")


# fit four primary glmm models and omnibus lrt tables ----
# note: format of data from pivot.confusion.metric.long
sens.analysis <- run.overall.confusion.analysis(
  long.data = sens.long.overall, 
  adjust.method = adj.method
  )

precision.analysis <- run.overall.confusion.analysis(
  long.data = precision.long.overall, 
  adjust.method = adj.method
  )

sens.bin.analysis <- run.bin.effects.confusion.analysis(
  long.data = sens.long.bins,
  adjust.method = adj.method
  )

precision.bin.analysis <- run.bin.effects.confusion.analysis(
  long.data = precision.long.bins,
  adjust.method = adj.method
  )

# build manuscript-facing plotting tables ----
plot.confusion.boxplots(sens.long, "ADMIX", "hapmap", "Sensitivity")
plot.confusion.boxplots(sens.long, "AFR", "hapmap", "Sensitivity")
plot.confusion.boxplots(sens.long, "ASIA", "hapmap", "Sensitivity")
plot.confusion.boxplots(sens.long, "EUR", "hapmap", "Sensitivity")

plot.confusion.boxplots(sens.long, "ADMIX", "decode", "Sensitivity")
plot.confusion.boxplots(sens.long, "AFR", "decode", "Sensitivity")
plot.confusion.boxplots(sens.long, "ASIA", "decode", "Sensitivity")
plot.confusion.boxplots(sens.long, "EUR", "decode", "Sensitivity")

plot.confusion.boxplots(sens.long, "ADMIX", "1Mb1cM", "Sensitivity")
plot.confusion.boxplots(sens.long, "AFR", "1Mb1cM", "Sensitivity")
plot.confusion.boxplots(sens.long, "ASIA", "1Mb1cM", "Sensitivity")
plot.confusion.boxplots(sens.long, "EUR", "1Mb1cM", "Sensitivity")

plot.confusion.boxplots(precision.long, "ADMIX", "hapmap", "Precision")
plot.confusion.boxplots(precision.long, "AFR", "hapmap", "Precision")
plot.confusion.boxplots(precision.long, "ASIA", "hapmap", "Precision")
plot.confusion.boxplots(precision.long, "EUR", "hapmap", "Precision")

plot.confusion.boxplots(precision.long, "ADMIX", "decode", "Precision")
plot.confusion.boxplots(precision.long, "AFR", "decode", "Precision")
plot.confusion.boxplots(precision.long, "ASIA", "decode", "Precision")
plot.confusion.boxplots(precision.long, "EUR", "decode", "Precision")

plot.confusion.boxplots(precision.long, "ADMIX", "1Mb1cM", "Precision")
plot.confusion.boxplots(precision.long, "AFR", "1Mb1cM", "Precision")
plot.confusion.boxplots(precision.long, "ASIA", "1Mb1cM", "Precision")
plot.confusion.boxplots(precision.long, "EUR", "1Mb1cM", "Precision")
