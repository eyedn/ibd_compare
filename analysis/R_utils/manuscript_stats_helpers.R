# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# manuscript_stats_helpers.R
# ______________________________________________________________________________


library(tidyverse)

# define shared manuscript model terms ----
manuscript.main.effects <- c("software.code", "filt", "imap", "pop")

manuscript.pairwise.interactions <- combn(
  manuscript.main.effects,
  2,
  FUN = function(x) paste(x, collapse = ":")
)

manuscript.full.fixed.terms <- c(
  manuscript.main.effects,
  manuscript.pairwise.interactions
)

# normalize compact directory-map strategy codes ----
normalize.manuscript.strategy.columns <- function(tbl) {
  if (!"software.code" %in% names(tbl)) {
    return(tbl)
  }

  if (!"imap" %in% names(tbl)) {
    tbl <- tbl %>%
      mutate(imap = NA_character_)
  }

  imap.code.map <- c(
    "h" = "hapmap",
    "d" = "decode",
    "f" = "1Mb1cM"
  )

  tbl %>%
    mutate(
      raw.software.code = as.character(software.code),
      strategy.match = str_match(raw.software.code, "^([hri])([hdf])(\\d+)$"),
      imap = if_else(
        !is.na(strategy.match[, 3]),
        unname(imap.code.map[strategy.match[, 3]]),
        as.character(.data[["imap"]])
      ),
      software.code = if_else(
        !is.na(strategy.match[, 2]),
        paste0(strategy.match[, 2], strategy.match[, 4]),
        raw.software.code
      )
    ) %>%
    select(-strategy.match)
}

# coerce model factors using the common manuscript columns ----
factor.manuscript.columns <- function(tbl) {
  factor.cols <- intersect(
    c("response", "bin", "software.code", "filt", "imap", "pop", "rep"),
    names(tbl)
  )

  tbl %>%
    mutate(across(all_of(factor.cols), factor))
}

# reshape confusion counts for sensitivity or precision models ----
build.manuscript.confusion.table <- function(all.stats.per.rep, response,
                                             bins) {
  stopifnot(response %in% c("sensitivity", "precision"))

  if (response == "sensitivity") {
    failure.prefix <- "fn"
    estimate.prefix <- "sensitivity"
  } else {
    failure.prefix <- "fp"
    estimate.prefix <- "ppv"
  }

  purrr::map_dfr(
    bins,
    function(bin) {
      suffix <- if (bin == "overall") "" else paste0(".", bin)
      tp.col <- paste0("tp", suffix)
      failure.col <- paste0(failure.prefix, suffix)
      estimate.col <- paste0(estimate.prefix, suffix)

      all.stats.per.rep %>%
        transmute(
          response = response,
          bin = bin,
          software.code = software.code,
          filt = filt,
          imap = imap,
          pop = pop,
          rep = rep,
          tp = .data[[tp.col]],
          failure = .data[[failure.col]],
          estimate = .data[[estimate.col]]
        )
    }
  ) %>%
    factor.manuscript.columns()
}

# build population-level IBD score error summaries ----
build.manuscript.ibd.score.table <- function(all.stats.per.rep,
                                             mse.epsilon) {
  all.stats.per.rep %>%
    transmute(
      response = "ibd.score",
      bin = "overall",
      software.code = software.code,
      filt = filt,
      imap = imap,
      pop = pop,
      rep = rep,
      true.ibd.score = true.ibd.score,
      inf.ibd.score = inf.ibd.score,
      bias = inf.ibd.score - true.ibd.score,
      mse = bias^2,
      log.mse = log(mse + mse.epsilon)
    ) %>%
    factor.manuscript.columns()
}

# summarize aligned upper-triangle kinship estimates by replicate strategy ----
build.manuscript.kinship.table <- function(true.kinship.tbl,
                                           inf.kinship.tbl,
                                           mse.epsilon,
                                           corr.clamp.epsilon) {
  true.tbl <- true.kinship.tbl %>%
    select(rep, pop, sam.pair, kinship) %>%
    rename(kin.true = kinship)

  inf.tbl <- inf.kinship.tbl %>%
    normalize.manuscript.strategy.columns() %>%
    select(rep, pop, sam.pair, software.code, filt, imap, kinship) %>%
    rename(kin.inf = kinship)

  inf.tbl %>%
    inner_join(true.tbl, by = c("rep", "pop", "sam.pair")) %>%
    filter(is.finite(kin.true), is.finite(kin.inf)) %>%
    group_by(rep, pop, software.code, filt, imap) %>%
    summarise(
      n.pairs = n(),
      bias = mean(kin.inf - kin.true),
      mse = mean((kin.inf - kin.true)^2),
      kinship.r = if_else(
        n() >= 2 && sd(kin.true) > 0 && sd(kin.inf) > 0,
        cor(kin.inf, kin.true),
        NA_real_
      ),
      .groups = "drop"
    ) %>%
    mutate(
      response = "kinship",
      bin = "overall",
      log.mse = log(mse + mse.epsilon),
      kinship.r.clamped = case_when(
        is.na(kinship.r) ~ NA_real_,
        is.finite(kinship.r) ~ pmin(
          pmax(kinship.r, -1 + corr.clamp.epsilon),
          1 - corr.clamp.epsilon
        ),
        TRUE ~ NA_real_
      ),
      kinship.z = atanh(kinship.r.clamped)
    ) %>%
    factor.manuscript.columns()
}

# convert raw KS p-values to manuscript model rows ----
build.manuscript.ks.table <- function(ks.tbl) {
  ks.tbl %>%
    normalize.manuscript.strategy.columns() %>%
    mutate(
      response = "ks.failure",
      bin = if ("bin" %in% names(.)) as.character(bin) else "overall",
      p.adj = p.adjust(p.value, method = "BH"),
      ks.failure = p.adj <= 0.05
    ) %>%
    select(
      response, bin, software.code, filt, imap, pop, rep,
      p.value, p.adj, ks.failure, everything()
    ) %>%
    factor.manuscript.columns()
}

# apply one global BH correction across all omnibus tests ----
add.global.bh.q.values <- function(omnibus.tbl) {
  omnibus.tbl %>%
    ungroup() %>%
    mutate(q.value = p.adjust(p.value, method = "BH"))
}

# choose pairwise terms retained in reduced final models ----
select.retained.interactions <- function(omnibus.tbl, q.threshold = 0.05) {
  omnibus.tbl %>%
    filter(str_detect(term, ":"), q.value < q.threshold) %>%
    pull(term) %>%
    unique()
}

# build model formulas from included fixed terms ----
build.manuscript.model.formula <- function(lhs, fixed.terms) {
  stats::as.formula(
    paste(
      lhs,
      "~",
      paste(fixed.terms, collapse = " + "),
      "+ (1 | rep)"
    )
  )
}

# build all omnibus reduced-term definitions ----
build.manuscript.omnibus.specs <- function() {
  main.specs <- purrr::map_dfr(
    manuscript.main.effects,
    function(term) {
      tibble(
        term = term,
        reduced.terms = list(
          manuscript.full.fixed.terms[
            !str_detect(
              manuscript.full.fixed.terms,
              paste0("(^|:)", term, "(:|$)")
            )
          ]
        )
      )
    }
  )

  interaction.specs <- purrr::map_dfr(
    manuscript.pairwise.interactions,
    function(term) {
      tibble(
        term = term,
        reduced.terms = list(setdiff(manuscript.full.fixed.terms, term))
      )
    }
  )

  bind_rows(main.specs, interaction.specs)
}

# build reduced final fixed terms from retained pairwise interactions ----
build.manuscript.reduced.fixed.terms <- function(retained.interactions) {
  c(
    manuscript.main.effects,
    intersect(manuscript.pairwise.interactions, retained.interactions)
  )
}
