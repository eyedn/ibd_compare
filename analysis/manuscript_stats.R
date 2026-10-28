# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# manuscript_stats.R
# ______________________________________________________________________________


library(tidyverse)
library(glue)
library(lme4)
library(emmeans)
library(furrr)
library(future)
library(progressr)

# source functions ----
analysis.dir <- if (dir.exists("analysis/R_utils")) {
  "analysis"
} else {
  "."
}

r.utils.dir <- "R_utils"
original.wd <- getwd()
setwd(analysis.dir)
on.exit(setwd(original.wd), add = TRUE)

source(file.path(r.utils.dir, "constants.R"))
source(file.path(r.utils.dir, "read.all.stats.per.rep.R"))
source(file.path(r.utils.dir, "read.pop.ibd.kinship.R"))
source(file.path(r.utils.dir, "read.pop.calls.R"))
source(file.path(r.utils.dir, "run.ks.R"))
source(file.path(r.utils.dir, "manuscript_stats_helpers.R"))

# configure runtime ----
cores <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")))
if (is.na(cores) || cores < 1) {
  cores <- 1
}
handlers("cli")

per.rep.path <- glue("{ibd.comp.dir}/true_vs_inference/all_stats.per_rep.tsv")
pop.info.dir <- glue("{sim.dir}/pop_info")

# local data readers ----
read.manuscript.unrels <- function(pop.info.dir) {
  tibble(
    file = list.files(
      pop.info.dir,
      full.names = TRUE,
      pattern = "^HapMapII_GRCh38_.*_unrel_samples\\.txt$"
    )
  ) %>%
    mutate(
      rep = as.integer(
        str_match(basename(file), "^HapMapII_GRCh38_(\\d+)_")[, 2]
      ),
      pop = str_match(
        basename(file),
        "^HapMapII_GRCh38_\\d+_([^_]+)_unrel_samples\\.txt$"
      )[, 2],
      ids = map(
        file,
        ~ readr::read_tsv(
          .x,
          col_names = FALSE,
          show_col_types = FALSE
        ) %>%
          transmute(ids = as.character(X1)) %>%
          pull(ids)
      )
    ) %>%
    select(rep, pop, ids)
}

read.manuscript.kinship.tables <- function(pops, dir.paths.map, cores) {
  true.kinship.tbl <- map_dfr(
    pops,
    ~ read.pop.ibd.kinship(.x, dir.paths.map, cores, which = "true")
  ) %>%
    select(-filt, -software.code) %>%
    distinct()

  inf.kinship.tbl <- map_dfr(
    pops,
    ~ read.pop.ibd.kinship(.x, dir.paths.map, cores, which = "inf")
  )

  list(
    true = true.kinship.tbl,
    inf = inf.kinship.tbl
  )
}

build.manuscript.ks.from.calls <- function(all.calls) {
  all.calls %>%
    normalize.manuscript.strategy.columns() %>%
    filter(software.code != "true", is.unrel == TRUE) %>%
    group_by(pop, rep, software.code, imap, filt) %>%
    summarise(inf = list(cM.len), .groups = "drop") %>%
    left_join(
      all.calls %>%
        filter(software.code == "true") %>%
        group_by(pop, rep) %>%
        summarise(tru = list(cM.len), .groups = "drop"),
      by = c("pop", "rep")
    ) %>%
    mutate(
      p.value = map2_dbl(tru, inf, run.ks),
      n.inf = map_int(inf, length),
      n.true = map_int(tru, length)
    ) %>%
    select(-tru, -inf) %>%
    build.manuscript.ks.table()
}

# local model helpers ----
fit.manuscript.model <- function(data, lhs, family = NULL,
                                 fixed.terms = manuscript.full.fixed.terms) {
  model.formula <- build.manuscript.model.formula(lhs, fixed.terms)

  if (is.null(family)) {
    lmer(model.formula, data = data, REML = FALSE)
  } else {
    glmer(model.formula, data = data, family = family)
  }
}

compare.manuscript.reduced.models <- function(data, lhs, family,
                                              full.model) {
  specs <- build.manuscript.omnibus.specs()

  map_dfr(
    seq_len(nrow(specs)),
    function(idx) {
      term <- specs$term[idx]
      reduced.terms <- specs$reduced.terms[[idx]]
      reduced.model <- fit.manuscript.model(
        data = data,
        lhs = lhs,
        family = family,
        fixed.terms = reduced.terms
      )
      anova.tbl <- as.data.frame(
        stats::anova(reduced.model, full.model, test = "Chisq")
      )

      tibble(
        term = term,
        reduced.formula = deparse(formula(reduced.model)),
        df.diff = anova.tbl$Df[2] - anova.tbl$Df[1],
        aic.diff = anova.tbl$AIC[1] - anova.tbl$AIC[2],
        lrt.statistic = anova.tbl$Chisq[2],
        p.value = anova.tbl$`Pr(>Chisq)`[2]
      )
    }
  )
}

extract.fixed.effects <- function(model) {
  coef.tbl <- as.data.frame(summary(model)$coefficients)
  coef.tbl$term <- rownames(coef.tbl)
  rownames(coef.tbl) <- NULL

  as_tibble(coef.tbl) %>%
    relocate(term)
}

extract.model.diagnostics <- function(model) {
  optinfo <- model@optinfo
  conv.messages <- optinfo$conv$lme4$messages

  tibble(
    is.singular = isSingular(model),
    optimizer = optinfo$optimizer,
    convergence.code = optinfo$conv$opt,
    convergence.message = paste(conv.messages, collapse = "; ")
  )
}

run.manuscript.posthoc <- function(model, retained.interactions,
                                   method.significant) {
  posthoc.specs <- tibble(
    term = c(
      "software.code:pop",
      "software.code:imap",
      "software.code:filt"
    ),
    formula = list(
      pairwise ~ software.code | pop,
      pairwise ~ software.code | imap,
      pairwise ~ software.code | filt
    )
  ) %>%
    filter(term %in% retained.interactions)

  interaction.res <- map_dfr(
    seq_len(nrow(posthoc.specs)),
    function(idx) {
      emm <- emmeans(model, posthoc.specs$formula[[idx]], adjust = "tukey")
      as_tibble(as.data.frame(emm$contrasts)) %>%
        mutate(posthoc.term = posthoc.specs$term[idx], .before = 1)
    }
  )

  if (nrow(interaction.res) > 0 || !method.significant) {
    return(interaction.res)
  }

  emm <- emmeans(model, pairwise ~ software.code, adjust = "tukey")
  as_tibble(as.data.frame(emm$contrasts)) %>%
    mutate(posthoc.term = "software.code", .before = 1)
}

# read and build manuscript model tables ----
all.stats.per.rep <- read.all.stats.per.rep(
  per.rep.path,
  stats.meta.names,
  param.map
)

sens.tbl <- build.manuscript.confusion.table(
  all.stats.per.rep = all.stats.per.rep,
  response = "sensitivity",
  bins = thr.levels
)
precision.tbl <- build.manuscript.confusion.table(
  all.stats.per.rep = all.stats.per.rep,
  response = "precision",
  bins = thr.levels
)
ibd.score.tbl <- build.manuscript.ibd.score.table(
  all.stats.per.rep = all.stats.per.rep,
  mse.epsilon = ibd.score.mse.epsilon
)

pops <- sort(unique(as.character(all.stats.per.rep$pop)))
kinship.raw <- read.manuscript.kinship.tables(pops, tvi.dir.paths.map, cores)
kinship.tbl <- build.manuscript.kinship.table(
  true.kinship.tbl = kinship.raw$true,
  inf.kinship.tbl = kinship.raw$inf,
  mse.epsilon = kinship.mse.epsilon,
  corr.clamp.epsilon = kinship.corr.clamp.epsilon
)

unrels.tbl <- read.manuscript.unrels(pop.info.dir)
all.calls <- map_dfr(
  pops,
  ~ read.pop.calls(.x, dir.paths.map, unrels.tbl, cores)
)
ks.tbl <- build.manuscript.ks.from.calls(all.calls)

manuscript.model.data <- list(
  sensitivity = sens.tbl,
  precision = precision.tbl,
  ibd.score = ibd.score.tbl,
  kinship = kinship.tbl,
  ks.failure = ks.tbl
)

# fit maximal models and omnibus likelihood-ratio tests ----
model.specs <- bind_rows(
  sens.tbl %>%
    distinct(response, bin) %>%
    mutate(lhs = "cbind(tp, failure)", family = list(binomial())),
  precision.tbl %>%
    distinct(response, bin) %>%
    mutate(lhs = "cbind(tp, failure)", family = list(binomial())),
  tibble(
    response = c(
      "ibd.score.bias", "ibd.score.log.mse",
      "kinship.bias", "kinship.log.mse", "kinship.correlation",
      "ks.failure"
    ),
    bin = "overall",
    lhs = c("bias", "log.mse", "bias", "log.mse", "kinship.z", "ks.failure"),
    family = list(NULL, NULL, NULL, NULL, NULL, binomial())
  )
) %>%
  mutate(model.id = paste(response, bin, sep = "__"))

model.specs <- model.specs %>%
  mutate(
    data = pmap(
      list(response, bin, lhs),
      function(response, bin, lhs) {
        response <- as.character(response)
        bin <- as.character(bin)

        if (str_starts(response, "sensitivity")) {
          sens.tbl %>% filter(bin == !!bin)
        } else if (str_starts(response, "precision")) {
          precision.tbl %>% filter(bin == !!bin)
        } else if (str_starts(response, "ibd.score")) {
          ibd.score.tbl
        } else if (str_starts(response, "kinship")) {
          kinship.tbl %>%
            filter(is.finite(.data[[lhs]]))
        } else {
          ks.tbl %>% filter(!is.na(ks.failure))
        }
      }
    )
  )

model.results <- pmap(
  list(model.specs$data, model.specs$lhs, model.specs$family),
  function(data, lhs, family) {
    full.model <- fit.manuscript.model(data, lhs, family)
    omnibus <- compare.manuscript.reduced.models(
      data = data,
      lhs = lhs,
      family = family,
      full.model = full.model
    )

    list(full.model = full.model, omnibus = omnibus)
  }
)
names(model.results) <- model.specs$model.id

manuscript.maximal.models <- map(model.results, "full.model")

manuscript.omnibus <- map2_dfr(
  model.results,
  model.specs$model.id,
  function(result, model.id) {
    model.specs %>%
      filter(model.id == !!model.id) %>%
      select(response, bin, model.id) %>%
      bind_cols(result$omnibus)
  }
) %>%
  add.global.bh.q.values()

# refit reduced final models with all main effects retained ----
reduced.model.specs <- model.specs %>%
  left_join(
    manuscript.omnibus %>%
      group_by(model.id) %>%
      summarise(
        retained.interactions = list(
          select.retained.interactions(
            tibble(term = term, q.value = q.value),
            q.threshold = 0.05
          )
        ),
        method.significant = any(
          term == "software.code" & q.value < 0.05,
          na.rm = TRUE
        ),
        .groups = "drop"
      ),
    by = "model.id"
  ) %>%
  mutate(
    reduced.fixed.terms = map(
      retained.interactions,
      build.manuscript.reduced.fixed.terms
    )
  )

manuscript.reduced.models <- pmap(
  list(
    reduced.model.specs$data,
    reduced.model.specs$lhs,
    reduced.model.specs$family,
    reduced.model.specs$reduced.fixed.terms
  ),
  fit.manuscript.model
)
names(manuscript.reduced.models) <- reduced.model.specs$model.id

manuscript.fixed.effects <- map2_dfr(
  manuscript.reduced.models,
  names(manuscript.reduced.models),
  function(model, model.id) {
    reduced.model.specs %>%
      filter(model.id == !!model.id) %>%
      select(response, bin, model.id) %>%
      bind_cols(extract.fixed.effects(model))
  }
)

manuscript.posthoc <- pmap_dfr(
  list(
    manuscript.reduced.models,
    reduced.model.specs$model.id,
    reduced.model.specs$retained.interactions,
    reduced.model.specs$method.significant
  ),
  function(model, model.id, retained.interactions, method.significant) {
    posthoc.tbl <- run.manuscript.posthoc(
      model = model,
      retained.interactions = retained.interactions,
      method.significant = method.significant
    )

    if (nrow(posthoc.tbl) == 0) {
      return(tibble())
    }

    reduced.model.specs %>%
      filter(model.id == !!model.id) %>%
      select(response, bin, model.id) %>%
      bind_cols(posthoc.tbl)
  }
)

manuscript.diagnostics <- map2_dfr(
  manuscript.reduced.models,
  names(manuscript.reduced.models),
  function(model, model.id) {
    reduced.model.specs %>%
      filter(model.id == !!model.id) %>%
      select(response, bin, model.id) %>%
      bind_cols(extract.model.diagnostics(model))
  }
)
