# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California 
# Department of Quantitative and Computational Biology 
# Mooney Lab 
# ___
# true_vs_inference_aggregate.R
# ______________________________________________________________________________

# libraries ----
library(tidyverse)
library(glue)

# source functions ----
source("job_scripts/R_utils/read.ibd.R")
source("job_scripts/R_utils/read.cov.res.R")
source("job_scripts/R_utils/compute.confusion.metrics.R")
source("job_scripts/R_utils/compute.confusion.metrics.windows.R")
source("job_scripts/R_utils/subset.unrels.R")
source("job_scripts/R_utils/compute.ibd.score.R")
source("job_scripts/R_utils/compute.kinship.coeffs.R")
source("job_scripts/R_utils/boot.mean.ci.R")
source("job_scripts/R_utils/infer.software.with.params.R")

# inputs ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7) {
  stop("Usage: true_vs_inference_aggregate.R MANIFEST.FILE TMAP SIM.CHR N.DIPLOIDS RES.OUT")
}
MANIFEST.FILE <- args[1]
TMAP <- args[2]
SIM.CHR <- args[3]
N.DIPLOIDS <- as.integer(args[4])
RES.OUT <- args[5]
FILT <- args[6]
POP <- args[7]

# constants ----
CHR1.100MB.LEN <- data.frame(
  "hapmap" = c(128.2627, 100000000),
  "flat"   = c(115.235, 100000000),
  "decode" = c(116.2396745, 100000000)
)
CHR18.LEN <- data.frame(
  "hapmap" = c(117.708923, 80257297),
  "flat"   = c(92.618155, 80373285),
  "decode" = c(105.783407, 80257297)
)

genome.seg.len.bp <- ifelse(SIM.CHR == "1", CHR1.100MB.LEN[2, TMAP], CHR18.LEN[2, TMAP])
genome.seg.len.cM <- ifelse(SIM.CHR == "1", CHR1.100MB.LEN[1, TMAP], CHR18.LEN[1, TMAP])

message("inputs:")
message("MANIFEST.FILE = ", MANIFEST.FILE)
message("TMAP = ", TMAP)
message("SIM.CHR = ", SIM.CHR)
message("N.DIPLOIDS = ", N.DIPLOIDS)
message("genome.seg.len.bp = ", genome.seg.len.bp)
message("genome.seg.len.cM = ", genome.seg.len.cM)
message("RES.OUT = ", RES.OUT)

# read manifest ----
manifest.df <- readr::read_tsv(MANIFEST.FILE, show_col_types = FALSE)
req.cols <- c("seed", "true_ibd", "inf_ibd", "coverage_out", "unrels")
if (!all(req.cols %in% names(manifest.df))) {
  stop("manifest must contain columns: ", paste(req.cols, collapse = ", "))
}

# parse metadata from manifest basename ----
# expected: software_tmap_imap_filt_pop.manifest.tsv
manifest.base <- basename(MANIFEST.FILE)
manifest.base <- sub("\\.manifest\\.tsv$", "", manifest.base)
parts <- strsplit(manifest.base, "_", fixed = TRUE)[[1]]

if (length(parts) < 5) {
  stop("manifest basename not parseable: ", basename(MANIFEST.FILE))
}

# infer software parameterization from directory name (software.dir.base)
software.base <- parts[1]
software.dir.base <- basename(dirname(MANIFEST.FILE))
software.with.params <- infer.software.with.params(software.base, software.dir.base)

message("meta data:")
message("software = ", parts[1])
message("software.with.params = ", software.with.params)
message("tmap = ", parts[2])
message("imap = ", parts[3])
message("filt = ", parts[4])
message("pop = ", parts[5])
meta.df <- tibble(
  "software" = parts[1],
  "software.with.params" = software.with.params,
  "tmap" = parts[2],
  "imap" = parts[3],
  "filt" = parts[4],
  "pop" = parts[5]
)

# output dir for per-rep kinship matrices ----
kin.dir <- file.path(dirname(MANIFEST.FILE), "ibd_kinship_mats")
dir.create(kin.dir, recursive = TRUE, showWarnings = FALSE)
message("kinship matrices will be saved to: ", kin.dir)

# per-replicate computations ----
n.reps <- nrow(manifest.df)
conf.list <- vector("list", n.reps)
ibd.list <- vector("list", n.reps)
K.true.list <- vector("list", n.reps)
K.inf.list  <- vector("list", n.reps)

for (i in seq_len(n.reps)) {
  
  seed <- manifest.df$seed[i]
  message("replicate ", i, "/", n.reps, " (seed=", seed, ")")
  
  true.file <- manifest.df$true_ibd[i]
  inf.file  <- manifest.df$inf_ibd[i]
  cov.file  <- manifest.df$coverage_out[i]
  unrel.file <- manifest.df$unrels[i]
  
  true.ibd <- read.ibd(true.file)
  inf.ibd  <- read.ibd(inf.file)
  cov.res  <- read.cov.res(cov.file)
  
  # confusion metrics (base + windowed) ----
  message("calculating confusion matrix statistics (overall)")
  conf.base <- compute.confusion.metrics(
    cov.res = cov.res,
    inf.ibd = inf.ibd,
    region.len.bp = genome.seg.len.bp,
    n.diploids = N.DIPLOIDS
  ) %>%
    tibble::as_tibble()

  message("calculating confusion matrix statistics (cM.len windows: 2-2.5, ..., 5+)")
  conf.win <- compute.confusion.metrics.windows(
    cov.res = cov.res,
    inf.ibd = inf.ibd,
    region.len.bp = genome.seg.len.bp,
    n.diploids = N.DIPLOIDS
  )

  conf.list[[i]] <- dplyr::bind_cols(conf.base, conf.win) %>%
    dplyr::mutate(seed = seed)

  # ibd scores (unrels)
  message("subsetting true ibd for unrels")
  true.unrel.res <- subset.unrels(unrel.file, true.ibd)
  message("subsetting inferred ibd for unrels")
  inf.unrel.res  <- subset.unrels(unrel.file, inf.ibd)
  
  message("calculating true ibd score")
  true.ibd.score <- compute.ibd.score(
    true.unrel.res$calls, true.unrel.res$n.unrels
    )
  message("calculating inferred ibd score")
  inf.ibd.score  <- compute.ibd.score(
    inf.unrel.res$calls,  inf.unrel.res$n.unrels
    )
  
  ibd.list[[i]] <- tibble(
    seed = seed,
    n.true.unrels = true.unrel.res$n.unrels,
    n.inf.unrels  = inf.unrel.res$n.unrels,
    true.ibd.score = true.ibd.score,
    inf.ibd.score  = inf.ibd.score,
    ibd.score.rel.true = inf.ibd.score / true.ibd.score,
    ibd.score.err = inf.ibd.score - true.ibd.score
  )
  
  # kinship matrices
  message("calculating true ibd-based kinship")
  K.true <- compute.kinship.coeffs(
    true.ibd, N.DIPLOIDS, genome.seg.len.cM
  )$kinship.matrix

  message("calculating inferred ibd-based kinship")
  K.inf <- compute.kinship.coeffs(
    inf.ibd,  N.DIPLOIDS, genome.seg.len.cM
  )$kinship.matrix

  # store in lists (unchanged behavior)
  K.true.list[[i]] <- K.true
  K.inf.list[[i]]  <- K.inf

  # write to disk for quick later access
  saveRDS(K.true, file = file.path(kin.dir, glue("true_{POP}_{seed}.{FILT}.rds")))
  saveRDS(K.inf,  file = file.path(kin.dir, glue("inf_{POP}_{seed}.{FILT}.rds")))
}

conf.df <- bind_rows(conf.list)
ibd.df  <- bind_rows(ibd.list)

# per-replicate output (no kinship) ----
message("writing per-replicate true-vs-inference summary (no kinship)")

per.rep.base <- conf.df %>%
  left_join(ibd.df, by = "seed") %>%
  mutate(rep = seed)

n.rows <- nrow(per.rep.base)

per.rep.df <- per.rep.base %>%
  bind_cols(meta.df[rep(1, n.rows), , drop = FALSE]) %>%
  select(
    software, software.with.params, tmap, imap, filt, pop,
    rep,
    # all IBD-related columns
    n.true.unrels, n.inf.unrels,
    true.ibd.score, inf.ibd.score,
    ibd.score.rel.true, ibd.score.err,
    # all confusion columns (whatever compute.confusion.metrics returned)
    all_of(setdiff(names(conf.df), "seed")),
    -seed
  )

PER_REP_OUT <- sub("\\.tsv$", ".per_rep.tsv", RES.OUT)
readr::write_tsv(per.rep.df, PER_REP_OUT, col_names = TRUE)
message("wrote per-replicate TSV: ", PER_REP_OUT)

# aggregate confusion metrics ----
message("aggregate confusion metrics")
conf.metrics <- setdiff(names(conf.df), "seed")

conf.summary <- map_dfc(conf.metrics, function(col) {
  ci <- boot.mean.ci(conf.df[[col]])
  tibble(
    !!paste0(col, ".mean") := ci["mean"],
    !!paste0(col, ".ci.lo") := ci["lo"],
    !!paste0(col, ".ci.hi") := ci["hi"]
  )
})

# aggregate ibd scores ----
message("aggregate ibd scores")
true.score.ci <- boot.mean.ci(ibd.df$true.ibd.score)
inf.score.ci  <- boot.mean.ci(ibd.df$inf.ibd.score)
rel.score.ci  <- boot.mean.ci(ibd.df$ibd.score.rel.true)

# bias/var/mse for (inf - true) ibd score across replicates
message("bias/var/mse for (inf - true) ibd score across replicates")
score.err <- ibd.df$ibd.score.err
ibd.score.bias <- mean(score.err, na.rm = TRUE)
ibd.score.var  <- var(score.err,  na.rm = TRUE)
ibd.score.mse  <- mean(score.err^2, na.rm = TRUE)

ibd.summary <- tibble(
  n.true.unrels.mean = mean(ibd.df$n.true.unrels),
  n.inf.unrels.mean  = mean(ibd.df$n.inf.unrels),
  
  true.ibd.score.mean = true.score.ci["mean"],
  true.ibd.score.ci.lo = true.score.ci["lo"],
  true.ibd.score.ci.hi = true.score.ci["hi"],
  
  inf.ibd.score.mean  = inf.score.ci["mean"],
  inf.ibd.score.ci.lo = inf.score.ci["lo"],
  inf.ibd.score.ci.hi = inf.score.ci["hi"],
  
  ibd.score.rel.true.mean  = rel.score.ci["mean"],
  ibd.score.rel.true.ci.lo = rel.score.ci["lo"],
  ibd.score.rel.true.ci.hi = rel.score.ci["hi"],
  
  ibd.score.bias = ibd.score.bias,
  ibd.score.var  = ibd.score.var,
  ibd.score.mse  = ibd.score.mse
)

# aggregate kinship: bias/var/mse across replicates, averaged over pairs ----
message("aggregate kinship: bias/var/mse across replicates")
idx <- upper.tri(K.true.list[[1]], diag = FALSE)

K.true.mat <- sapply(K.true.list, function(K) K[idx])  # [n.pairs, n.reps]
K.inf.mat  <- sapply(K.inf.list,  function(K) K[idx])

K.inf.bar  <- rowMeans(K.inf.mat)
K.true.bar <- rowMeans(K.true.mat)

bias.ij <- K.inf.bar - K.true.bar
var.ij  <- apply(K.inf.mat, 1, var)
mse.ij  <- rowMeans((K.inf.mat - K.true.mat)^2)

kinship.summary <- tibble(
  kinship.bias = mean(bias.ij),
  kinship.var  = mean(var.ij),
  kinship.mse  = mean(mse.ij)
)

# final one-row output ----
res.df <- bind_cols(
  meta.df,
  tibble(n.reps = n.reps),
  conf.summary,
  ibd.summary,
  kinship.summary
)

readr::write_tsv(res.df, RES.OUT, col_names = TRUE)
message("wrote: ", RES.OUT)
