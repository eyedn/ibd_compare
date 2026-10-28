# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# zero_inferred_ibd.R
# ______________________________________________________________________________

# Run from the repository root: Rscript tests/zero_inferred_ibd.R

suppressPackageStartupMessages(library(tidyverse))
source("job_scripts/R_utils/read.ibd.R")
source("job_scripts/R_utils/compute.ibd.score.R")
source("job_scripts/R_utils/subset.unrels.R")
source("job_scripts/R_utils/compute.kinship.coeffs.R")

fixture.dir <- tempfile("zero-inferred-ibd-", tmpdir = "/tmp")
dir.create(fixture.dir)
failures <- character()

check <- function(label, code) {
  tryCatch({
    force(code)
    message("PASS: ", label)
  }, error = function(err) {
    failures <<- c(failures, label)
    message("FAIL: ", label, ": ", conditionMessage(err))
  })
}

expect.failure <- function(path, pattern) {
  # Invalid inputs intentionally trigger readr diagnostics before our error.
  err <- suppressWarnings(tryCatch(read.ibd(path), error = identity))
  stopifnot(inherits(err, "error"))
  stopifnot(grepl(path, conditionMessage(err), fixed = TRUE))
  stopifnot(grepl(pattern, conditionMessage(err), ignore.case = TRUE))
}

# A single 5 cM segment in a three-sample population.
calls <- tibble(
  chr = "1", bp.str = 100L, bp.end = 1100L, bp.len = 1000L,
  cM.str = 1, cM.end = 6, cM.len = 5,
  sam1 = "TEST_1", hap1 = 1L, sam2 = "TEST_2", hap2 = 2L
)
empty.file <- file.path(fixture.dir, "empty.bed.gz")
full.file <- file.path(fixture.dir, "full.bed.gz")
write_tsv(calls[0, ], empty.file)
write_tsv(calls, full.file)
unrel.file <- file.path(fixture.dir, "unrels.txt")
writeLines(paste0("TEST_", 1:3), unrel.file)

check("header-only gzip has standard types and score zero", {
  empty <- read.ibd(empty.file)
  stopifnot(nrow(empty) == 0L, identical(names(empty), names(calls)))
  stopifnot(identical(vapply(empty, typeof, character(1)),
                      vapply(calls, typeof, character(1))))
  stopifnot(identical(compute.ibd.score(empty, 3), 0))
})

check("unrelated-pair filtering can remove every segment", {
  excluded <- mutate(calls, sam2 = "TEST_4")
  unrelated <- subset.unrels(unrel.file, excluded)
  stopifnot(nrow(unrelated$calls) == 0L)
  stopifnot(identical(compute.ibd.score(unrelated$calls, 3), 0))
})

check("length filtering can remove every segment", {
  short <- mutate(calls, cM.len = 2)
  stopifnot(identical(compute.ibd.score(short, 3), 0))
})

check("nonempty score, inclusive limits, and kinship stay unchanged", {
  full <- read.ibd(full.file)
  stopifnot(compute.ibd.score(full, 3) == 1000 / 12)
  limits <- bind_rows(mutate(full, cM.len = 3),
                      mutate(full, cM.len = 20))
  stopifnot(compute.ibd.score(limits, 3) == 2000 / 12)
  expected <- diag(0.5, 3)
  dimnames(expected) <- rep(list(paste0("TEST_", 1:3)), 2)
  expected[1, 2] <- expected[2, 1] <- 5 / (4 * 128.2627)
  actual <- compute.kinship.coeffs(
    full, 3, 128.2627, sample.ids = paste0("TEST_", 1:3)
  )$kinship.matrix
  stopifnot(isTRUE(all.equal(actual, expected)))
})

check("missing columns are rejected with filename", {
  path <- file.path(fixture.dir, "missing-columns.bed.gz")
  write_tsv(select(calls, -hap2), path)
  expect.failure(path, "missing.*hap2")
})

for (column in c("bp.str", "bp.end", "bp.len", "hap1", "hap2",
                 "cM.str", "cM.end", "cM.len")) {
  check(paste("malformed numeric field is rejected:", column), {
    path <- file.path(fixture.dir, paste0("malformed-", column, ".bed.gz"))
    malformed <- calls
    malformed[[column]] <- "bad-number"
    write_tsv(malformed, path)
    expect.failure(path, "pars")
  })
}

check("headerless input is rejected with filename", {
  path <- file.path(fixture.dir, "headerless.bed.gz")
  write_tsv(calls, path, col_names = FALSE)
  expect.failure(path, "missing|header|column")
})

check("completely empty input is rejected with filename", {
  path <- file.path(fixture.dir, "no-header.bed")
  file.create(path)
  expect.failure(path, "missing|header|column")
})

check("missing input file remains an error", {
  expect.failure(file.path(fixture.dir, "nonexistent.bed.gz"),
                 "exist|open|file")
})

# Freeze the public output schemas independently of the script under test.
meta.names <- c("software", "software.with.params", "tmap", "imap", "filt",
                "pop")
metric.names <- c("tp", "fn", "fp", "tn", "sensitivity", "specificity",
                  "ppv", "npv", "fpr", "accuracy")
conf.names <- c(metric.names, unlist(lapply(
  c(seq(2, 9.5, by = 0.5), 10), function(lo) paste0(metric.names, ".", lo)
)))
per.rep.names <- c(
  meta.names, "rep", "n.true.unrels", "n.inf.unrels", "true.ibd.score",
  "inf.ibd.score", "ibd.score.rel.true", "ibd.score.err", conf.names
)
agg.names <- c(
  meta.names, "n.reps",
  unlist(lapply(conf.names, function(col) {
    paste0(col, c(".mean", ".ci.lo", ".ci.hi"))
  })),
  "n.true.unrels.mean", "n.inf.unrels.mean",
  paste0("true.ibd.score", c(".mean", ".ci.lo", ".ci.hi")),
  paste0("inf.ibd.score", c(".mean", ".ci.lo", ".ci.hi")),
  paste0("ibd.score.rel.true", c(".mean", ".ci.lo", ".ci.hi")),
  "ibd.score.bias", "ibd.score.var", "ibd.score.mse",
  "kinship.bias", "kinship.var", "kinship.mse"
)

# Exercise both manifest names with mixed, all-empty, and nonempty controls.
for (software in c("rapid", "hapibd")) {
  for (scenario in c("mixed", "all-empty", "nonempty")) {
    check(paste("aggregation", software, scenario), {
      run.dir <- file.path(fixture.dir, scenario, software)
      dir.create(run.dir, recursive = TRUE)
      is.empty <- switch(
        scenario, mixed = c(TRUE, FALSE), `all-empty` = c(TRUE, TRUE),
        nonempty = c(FALSE, FALSE)
      )
      cov.files <- file.path(run.dir, paste0("coverage-", 1:2, ".tsv"))
      for (i in 1:2) {
        coverage <- mutate(
          calls, num.infer.ibd = as.integer(!is.empty[i]),
          bp.cov = if (is.empty[i]) 0L else 1000L,
          bp.perc.cov = if (is.empty[i]) 0 else 1
        )
        write_tsv(coverage, cov.files[i], col_names = FALSE)
      }
      manifest <- file.path(
        run.dir, paste0(software, "_hapmap_hapmap_filt_TEST.manifest.tsv")
      )
      write_tsv(tibble(
        rep = 1:2, true_ibd = full.file,
        inf_ibd = ifelse(is.empty, empty.file, full.file),
        coverage_out = cov.files, unrels = unrel.file
      ), manifest)
      output <- file.path(run.dir, "results.agg.tsv")
      log.file <- file.path(run.dir, "aggregation.log")
      status <- system2(file.path(R.home("bin"), "Rscript"), c(
        "job_scripts/true_vs_inference_aggregate.R", shQuote(manifest),
        "hapmap", "1", "3", shQuote(output), "filt", "TEST"
      ), stdout = log.file, stderr = log.file)
      if (status != 0) stop(paste(readLines(log.file), collapse = "\n"))
      per.rep <- read_tsv(sub(".tsv$", ".per_rep.tsv", output),
                          show_col_types = FALSE)
      agg <- read_tsv(output, show_col_types = FALSE)
      stopifnot(nrow(per.rep) == 2, nrow(agg) == 1, agg$n.reps == 2)
      stopifnot(identical(names(per.rep), per.rep.names))
      stopifnot(identical(names(agg), agg.names))
      stopifnot(identical(per.rep$rep, as.double(1:2)))
      stopifnot(all(per.rep$true.ibd.score == 1000 / 12))
      scores <- ifelse(is.empty, 0, 1000 / 12)
      stopifnot(isTRUE(all.equal(per.rep$inf.ibd.score, scores)))
      stopifnot(isTRUE(all.equal(per.rep$ibd.score.err, scores - 1000 / 12)))
      stopifnot(all(per.rep$ibd.score.rel.true == ifelse(is.empty, 0, 1)))
      stopifnot(all(per.rep$tp == ifelse(is.empty, 0, 1000)))
      stopifnot(all(per.rep$fn == ifelse(is.empty, 1000, 0)))
      stopifnot(all(per.rep$fp == 0), all(per.rep$tn == 1200000000 - 1000))
      stopifnot(all(is.na(per.rep$ppv[is.empty])))
      stopifnot(all(per.rep$ppv[!is.empty] == 1))
      if (all(is.empty)) stopifnot(is.na(agg$ppv.mean))
      stopifnot(isTRUE(all.equal(agg$inf.ibd.score.mean, mean(scores))))
      stopifnot(isTRUE(all.equal(agg$ibd.score.bias,
                                mean(scores - 1000 / 12))))
      kinship.err <- ifelse(is.empty, -5 / (4 * 128.2627), 0)
      stopifnot(isTRUE(all.equal(agg$kinship.bias, mean(kinship.err) / 3)))
      stopifnot(isTRUE(all.equal(agg$kinship.var, var(kinship.err) / 3)))
      stopifnot(isTRUE(all.equal(agg$kinship.mse,
                                mean(kinship.err^2) / 3)))
      log.lines <- readLines(log.file)
      for (i in 1:2) {
        expected <- diag(0.5, 3)
        dimnames(expected) <- rep(list(paste0("TEST_", 1:3)), 2)
        expected[1, 2] <- expected[2, 1] <- 5 / (4 * 128.2627)
        for (kind in c("true", "inf")) {
          path <- file.path(run.dir, "ibd_kinship_mats",
                            paste0(kind, "_TEST_", i, ".filt.rds"))
          target <- expected
          if (kind == "inf" && is.empty[i]) {
            target[row(target) != col(target)] <- 0
          }
          stopifnot(isTRUE(all.equal(readRDS(path), target)))
        }
        if (is.empty[i]) {
          log.match <- log.lines[grepl("zero.*segments", log.lines)]
          stopifnot(any(grepl(paste0("rep=", i), log.match, fixed = TRUE) &
                          grepl(empty.file, log.match, fixed = TRUE)))
        }
      }
    })
  }
}

if (length(failures)) {
  stop(length(failures), " checks failed; fixtures retained at ", fixture.dir)
}
unlink(fixture.dir, recursive = TRUE)
message("All zero-inferred-IBD regression checks passed (local fixtures only).")
