# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# calc.3T.Ne.R
# ______________________________________________________________________________


library(tidyverse)

# functions ----
# compute the founder effective pop size from a founder census size.
# `founder.count` is the total founder count and `female.proportion` is the
# fraction of founders that are female.
calc.founder.ne <- function(founder.count, female.proportion) {
  female.count <- founder.count * female.proportion
  male.count <- founder.count - female.count

  if (female.count <= 0 || male.count <= 0) {
    return(0)
  }

  4 * female.count * male.count / (female.count + male.count)
}

# validate that a given value is a single non-negative whole generation.
# `value` is the candidate input and `name` is used in the error message.
validate.whole.generation.count <- function(value, name) {
  if (
    length(value) != 1 ||
      is.na(value) ||
      !is.finite(value) ||
      value != as.integer(value) ||
      value < 0
  ) {
    stop(name, " must be a non-negative whole generation count", call. = FALSE)
  }

  as.integer(value)
}

# build a generation sequence that can be empty when the requested length is 0.
# `start` is the first generation in the block and `len` is the block length.
build.generation.block <- function(start, len) {
  if (len == 0) {
    return(integer(0))
  }

  seq.int(start, start + len - 1L)
}

# return an explicit argument when provided, otherwise pull the default.
# `value` is the user-supplied override, `name` is the default binding name,
# and `defaults.env` is the environment that stores the fallback value.
resolve.default.param <- function(value, name, defaults.env) {
  if (!is.null(value)) {
    return(value)
  }

  get(name, envir = defaults.env)
}

# build the generation layout for the three epochs.
# `epoch.*.len` are the epoch lengths, `gen.time` and `present.year` define the
#  timeline, and `start.year` sets the first generation boundary.
calc.epoch.config <- function(
    epoch.1.len,
    epoch.2.len,
    epoch.3.len,
    gen.time,
    present.year,
    start.year = 1600
) {
  epoch.1.len <- validate.whole.generation.count(
    epoch.1.len,
    "epoch.1.len"
  )
  epoch.2.len <- validate.whole.generation.count(
    epoch.2.len,
    "epoch.2.len"
  )
  epoch.3.len <- validate.whole.generation.count(
    epoch.3.len,
    "epoch.3.len"
  )

  year.span <- present.year - start.year
  total.generations <- year.span / gen.time

  if (total.generations != as.integer(total.generations)) {
    stop(
      "present.year - 1600 must be divisible by gen.time",
      call. = FALSE
    )
  }

  total.generations <- as.integer(total.generations)
  total.epoch.len <- epoch.1.len + epoch.2.len + epoch.3.len

  if (total.epoch.len != total.generations) {
    stop(
      "total epoch length must equal (present.year - 1600) / gen.time",
      call. = FALSE
    )
  }

  epoch.1.generations <- build.generation.block(1L, epoch.1.len)
  epoch.2.generations <- build.generation.block(epoch.1.len + 1L, epoch.2.len)
  epoch.3.first.generation <- if (epoch.3.len == 0) {
    NA_integer_
  } else {
    epoch.1.len + epoch.2.len + 1L
  }
  epoch.3.generations <- build.generation.block(
    epoch.1.len + epoch.2.len + 1L,
    epoch.3.len
  )

  epoch.3.anchor.generation <- NA_integer_

  if (epoch.3.len > 0) {
    if (epoch.2.len > 0) {
      epoch.3.anchor.generation <- max(epoch.2.generations)
    } else if (epoch.1.len > 0) {
      epoch.3.anchor.generation <- max(epoch.1.generations)
    } else {
      stop(
        "epoch 3 requires a prior non-empty epoch to provide an anchor",
        call. = FALSE
      )
    }
  }

  list(
    epoch.1.len = epoch.1.len,
    epoch.2.len = epoch.2.len,
    epoch.3.len = epoch.3.len,
    total.generations = total.generations,
    final.generation = total.generations,
    epoch.1.generations = epoch.1.generations,
    epoch.1.sex.generations = epoch.1.generations[epoch.1.generations >= 2],
    epoch.2.generations = epoch.2.generations,
    epoch.2.parental.generations = epoch.2.generations[
      epoch.2.generations >= 2
    ],
    epoch.3.generations = epoch.3.generations,
    epoch.3.first.generation = epoch.3.first.generation,
    epoch.3.anchor.generation = epoch.3.anchor.generation
  )
}

# create the African and European ancestry proportions by generation.
# `epoch.config` controls the generation grid, `mooney.cg` the proportion of
# non-European ancestry that is African in epoch-1, and the remaining arguments
# define the fixed mixture proportions in each following epoch.
build.mixing.props <- function(
    epoch.config,
    mooney.cg,
    epoch.1.eur.prop,
    epoch.2.afr.prop,
    epoch.2.eur.prop,
    epoch.3.afr.prop,
    epoch.3.eur.prop
) {
  if (epoch.config$epoch.1.len > length(mooney.cg)) {
    stop("epoch.1.len cannot exceed length(mooney.cg)", call. = FALSE)
  }

  s.eur <- c(
    rep(epoch.1.eur.prop, epoch.config$epoch.1.len),
    rep(epoch.2.eur.prop, epoch.config$epoch.2.len),
    rep(epoch.3.eur.prop, epoch.config$epoch.3.len)
  )
  s.afr <- c(
    (1 - epoch.1.eur.prop) * mooney.cg[seq_len(epoch.config$epoch.1.len)],
    rep(epoch.2.afr.prop, epoch.config$epoch.2.len),
    rep(epoch.3.afr.prop, epoch.config$epoch.3.len)
  )

  if (
    length(s.eur) != epoch.config$total.generations ||
      length(s.afr) != epoch.config$total.generations
  ) {
    stop("AA.mixing.props vectors must not recycle", call. = FALSE)
  }

  tibble(
    g = seq_len(epoch.config$total.generations),
    s.afr = s.afr,
    s.eur = s.eur,
    h = pmax(1 - s.afr - s.eur, 0)
  )
}

# intersect generation windows with historical periods and measure coverage.
# `window.tbl` is the generation map and `hacker.tbl` is the historical  period
# table with start and end years.
calc.period.overlaps <- function(window.tbl, hacker.tbl) {
  merge(window.tbl, hacker.tbl, by = NULL) %>%
    as_tibble() %>%
    mutate(
      period.length = end.year - start.year,
      overlap.start.year = pmax(window.start.year, start.year),
      overlap.end.year = pmin(window.end.year, end.year),
      overlap.years = pmax(0, overlap.end.year - overlap.start.year),
      overlap.weight = overlap.years / period.length
    ) %>%
    filter(overlap.years > 0) %>%
    arrange(generation, start.year)
}

# estimate population size at a year by interpolating inside a period.
# `year` is the target year to query and `hacker.tbl` contains the period
# boundaries plus beginning and ending population sizes.
calc.boundary.population <- function(year, hacker.tbl) {
  if (is.na(year)) {
    return(NA_real_)
  }

  if (year < min(hacker.tbl$start.year) || year > max(hacker.tbl$end.year)) {
    return(NA_real_)
  }

  period.idx <- which(
    hacker.tbl$start.year <= year & hacker.tbl$end.year >= year
  )
  period.idx <- period.idx[length(period.idx)]
  period.row <- hacker.tbl[period.idx, ]
  period.length <- period.row$end.year - period.row$start.year
  period.frac <- (year - period.row$start.year) / period.length

  period.row$pop.beginning +
    period.frac * (period.row$pop.end - period.row$pop.beginning)
}

# add a readable epoch label to each generation row.
# `AA.ne` is the analysis table and `epoch.config` supplies the generation
# ranges used to assign the three epoch labels.
make.epoch.labels <- function(AA.ne, epoch.config) {
  epoch.1.label <- if (length(epoch.config$epoch.1.generations) > 0) {
    paste0(
      "Epoch 1: G",
      min(epoch.config$epoch.1.generations),
      "-G",
      max(epoch.config$epoch.1.generations),
      ", Import-influenced sex-ratio"
    )
  } else {
    NA_character_
  }
  epoch.2.label <- if (length(epoch.config$epoch.2.generations) > 0) {
    paste0(
      "Epoch 2: G",
      min(epoch.config$epoch.2.generations),
      "-G",
      max(epoch.config$epoch.2.generations),
      ", WF assumptions"
    )
  } else {
    NA_character_
  }
  epoch.3.label <- if (length(epoch.config$epoch.3.generations) > 0) {
    paste0(
      "Epoch 3: G",
      min(epoch.config$epoch.3.generations),
      "-G",
      max(epoch.config$epoch.3.generations),
      ", modern exponential growth"
    )
  } else {
    NA_character_
  }

  AA.ne %>%
    mutate(
      epoch = case_when(
        generation %in% epoch.config$epoch.1.generations ~ epoch.1.label,
        generation %in% epoch.config$epoch.2.generations ~ epoch.2.label,
        generation %in% epoch.config$epoch.3.generations ~ epoch.3.label
      )
    )
}

# check that the derived tables are internally consistent before returning.
# `AA.ne`, `AA.overlap`, and `AA.mixing.props` are the computed outputs, while
# the remaining arguments provide the timeline and boundary checks.
validate.analysis.outputs <- function(
    AA.ne,
    AA.overlap,
    AA.mixing.props,
    epoch.config,
    hacker.start.year,
    gen.time,
    present.year
) {
  stopifnot(
    identical(
      sort(AA.ne$generation),
      seq_len(epoch.config$total.generations)
    ),
    nrow(AA.mixing.props) == epoch.config$total.generations,
    identical(AA.mixing.props$g, seq_len(epoch.config$total.generations)),
    all(AA.ne$generation.window.years == gen.time),
    AA.ne$aggregation.start.year[1] == hacker.start.year,
    AA.ne$generation.end.year[epoch.config$final.generation] == present.year,
    all(AA.overlap$overlap.weight > 0),
    all(AA.overlap$overlap.weight <= 1),
    all(!is.na(AA.ne$female.total[epoch.config$epoch.1.sex.generations])),
    all(!is.na(AA.ne$male.total[epoch.config$epoch.1.sex.generations])),
    all(
      !is.na(
        AA.ne$parental.population[epoch.config$epoch.2.parental.generations]
      )
    ),
    all(!is.na(AA.ne$k[epoch.config$epoch.2.generations])),
    all(
      grepl(
        "^Epoch 3:",
        AA.ne$epoch[epoch.config$epoch.3.generations]
      )
    ),
    all(is.finite(AA.ne$ne)),
    all(AA.ne$ne > 0)
  )
}

# compute the full AA Ne analysis and return every intermediate table.
# the scalar arguments control the demographic timeline, founders, and epoch
# settings, while `mooney.cg` and `hacker.tbl` provide the input data.
calc.aa.ne.results <- function(
    gen.time = NULL,
    hacker.start.year = NULL,
    hacker.end.year = NULL,
    present.year = NULL,
    epoch.1.len = NULL,
    epoch.2.len = NULL,
    epoch.3.len = NULL,
    afr.founder.count = NULL,
    afr.female.proportion = NULL,
    eur.founder.count = NULL,
    eur.female.proportion = NULL,
    epoch.3.growth.rate = NULL,
    mooney.cg = NULL,
    epoch.1.eur.prop = NULL,
    epoch.2.afr.prop = NULL,
    epoch.2.eur.prop = NULL,
    epoch.3.afr.prop = NULL,
    epoch.3.eur.prop = NULL,
    hacker.tbl = NULL
) {
  defaults.env <- parent.env(environment())
  gen.time <- resolve.default.param(gen.time, "gen.time", defaults.env)
  hacker.start.year <- resolve.default.param(
    hacker.start.year,
    "hacker.start.year",
    defaults.env
  )
  hacker.end.year <- resolve.default.param(
    hacker.end.year,
    "hacker.end.year",
    defaults.env
  )
  present.year <- resolve.default.param(
    present.year,
    "present.year",
    defaults.env
  )
  epoch.1.len <- resolve.default.param(
    epoch.1.len,
    "epoch.1.len",
    defaults.env
  )
  epoch.2.len <- resolve.default.param(
    epoch.2.len,
    "epoch.2.len",
    defaults.env
  )
  epoch.3.len <- resolve.default.param(
    epoch.3.len,
    "epoch.3.len",
    defaults.env
  )
  afr.founder.count <- resolve.default.param(
    afr.founder.count,
    "afr.founder.count",
    defaults.env
  )
  afr.female.proportion <- resolve.default.param(
    afr.female.proportion,
    "afr.female.proportion",
    defaults.env
  )
  eur.founder.count <- resolve.default.param(
    eur.founder.count,
    "eur.founder.count",
    defaults.env
  )
  eur.female.proportion <- resolve.default.param(
    eur.female.proportion,
    "eur.female.proportion",
    defaults.env
  )
  epoch.3.growth.rate <- resolve.default.param(
    epoch.3.growth.rate,
    "epoch.3.growth.rate",
    defaults.env
  )
  mooney.cg <- resolve.default.param(mooney.cg, "mooney.cg", defaults.env)
  epoch.1.eur.prop <- resolve.default.param(
    epoch.1.eur.prop,
    "epoch.1.eur.prop",
    defaults.env
  )
  epoch.2.afr.prop <- resolve.default.param(
    epoch.2.afr.prop,
    "epoch.2.afr.prop",
    defaults.env
  )
  epoch.2.eur.prop <- resolve.default.param(
    epoch.2.eur.prop,
    "epoch.2.eur.prop",
    defaults.env
  )
  epoch.3.afr.prop <- resolve.default.param(
    epoch.3.afr.prop,
    "epoch.3.afr.prop",
    defaults.env
  )
  epoch.3.eur.prop <- resolve.default.param(
    epoch.3.eur.prop,
    "epoch.3.eur.prop",
    defaults.env
  )
  hacker.tbl <- resolve.default.param(hacker.tbl, "hacker.tbl", defaults.env)

  epoch.config <- calc.epoch.config(
    epoch.1.len = epoch.1.len,
    epoch.2.len = epoch.2.len,
    epoch.3.len = epoch.3.len,
    gen.time = gen.time,
    present.year = present.year
  )

  AA.mixing.props <- build.mixing.props(
    epoch.config = epoch.config,
    mooney.cg = mooney.cg,
    epoch.1.eur.prop = epoch.1.eur.prop,
    epoch.2.afr.prop = epoch.2.afr.prop,
    epoch.2.eur.prop = epoch.2.eur.prop,
    epoch.3.afr.prop = epoch.3.afr.prop,
    epoch.3.eur.prop = epoch.3.eur.prop
  )

  AA.hacker.tbl <- hacker.tbl %>%
    mutate(
      start.year = as.numeric(sub("-.*", "", period)),
      end.year = as.numeric(sub(".*-", "", period))
    ) %>%
    select(
      period,
      start.year,
      end.year,
      person.years.lived,
      assumed.birth.rate,
      imported,
      births,
      pop.beginning,
      pop.end
    )

  AA.generation.map <- tibble(
    generation = seq_len(epoch.config$total.generations),
    generation.start.year = seq(
      1600,
      present.year - gen.time,
      by = gen.time
    ),
    generation.end.year = generation.start.year + gen.time,
    generation.window.years = generation.end.year - generation.start.year,
    window.start.year = generation.start.year,
    window.end.year = generation.end.year
  )

  AA.overlap <- calc.period.overlaps(
    AA.generation.map,
    AA.hacker.tbl
  ) %>%
    select(
      generation,
      generation.start.year,
      generation.end.year,
      period,
      start.year,
      end.year,
      overlap.start.year,
      overlap.end.year,
      overlap.years,
      overlap.weight,
      imported,
      births,
      assumed.birth.rate,
      person.years.lived
    )

  AA.agg <- AA.overlap %>%
    group_by(generation, generation.start.year, generation.end.year) %>%
    summarise(
      periods = paste(period, collapse = ", "),
      imported = sum(imported * overlap.weight, na.rm = TRUE),
      births = sum(births * overlap.weight, na.rm = TRUE),
      birth.rate = weighted.mean(
        assumed.birth.rate,
        person.years.lived * overlap.weight,
        na.rm = TRUE
      ),
      aggregation.start.year = min(overlap.start.year),
      aggregation.end.year = max(overlap.end.year),
      hacker.window.years = sum(overlap.years),
      .groups = "drop"
    )

  AA.agg <- AA.generation.map %>%
    select(
      generation,
      generation.start.year,
      generation.end.year,
      generation.window.years
    ) %>%
    left_join(
      AA.agg,
      by = c("generation", "generation.start.year", "generation.end.year")
    ) %>%
    mutate(
      population.boundary.year = if_else(
        generation.start.year < hacker.end.year,
        pmin(generation.end.year, hacker.end.year),
        NA_real_
      ),
      pop.end = map_dbl(
        population.boundary.year,
        calc.boundary.population,
        hacker.tbl = AA.hacker.tbl
      ),
      is.truncated = !is.na(hacker.window.years) &
        hacker.window.years < generation.window.years
    )

  AA.ne <- AA.agg %>%
    arrange(generation) %>%
    mutate(
      population.at.generation.start = map_dbl(
        generation.start.year,
        calc.boundary.population,
        hacker.tbl = AA.hacker.tbl
      ),
      prev.births = dplyr::lag(births),
      previous.birth.periods = NA_character_,
      births.in.previous.25.years = NA_real_,
      female.imported = if_else(
        generation %in% epoch.config$epoch.1.sex.generations,
        ceiling(imported / 3),
        NA_real_
      ),
      male.imported = if_else(
        generation %in% epoch.config$epoch.1.sex.generations,
        floor(2 * imported / 3),
        NA_real_
      ),
      female.potential = if_else(
        generation %in% epoch.config$epoch.1.sex.generations,
        0.5 * prev.births,
        NA_real_
      ),
      male.potential = if_else(
        generation %in% epoch.config$epoch.1.sex.generations,
        0.5 * prev.births,
        NA_real_
      ),
      female.total = female.imported + female.potential,
      male.total = male.imported + male.potential,
      adult.population.age25.proxy = NA_real_,
      parental.population = if_else(
        generation %in% epoch.config$epoch.2.generations,
        prev.births,
        NA_real_
      ),
      parental.population.method = if_else(
        generation %in% epoch.config$epoch.2.generations,
        "previous-generation births",
        NA_character_
      ),
      k = if_else(
        generation %in% epoch.config$epoch.2.generations,
        gen.time * birth.rate / 1000,
        NA_real_
      )
    )

  ne.g1 <- calc.founder.ne(
    afr.founder.count,
    afr.female.proportion
  ) +
    calc.founder.ne(
      eur.founder.count,
      eur.female.proportion
    )

  AA.ne <- AA.ne %>%
    mutate(
      ne = case_when(
        generation == 1 ~ ne.g1,
        generation %in% epoch.config$epoch.1.sex.generations ~
          4 * female.total * male.total / (female.total + male.total),
        generation %in% epoch.config$epoch.2.generations ~
          (k * parental.population - 1) / k,
        TRUE ~ NA_real_
      )
    )

  epoch.3.anchor.ne <- if (length(epoch.config$epoch.3.generations) == 0) {
    NA_real_
  } else {
    AA.ne %>%
      filter(generation == epoch.config$epoch.3.anchor.generation) %>%
      pull(ne)
  }

  AA.ne <- AA.ne %>%
    mutate(
      ne = if_else(
        generation %in% epoch.config$epoch.3.generations,
        epoch.3.anchor.ne *
          exp(epoch.3.growth.rate *
            (generation - epoch.config$epoch.3.anchor.generation)),
        ne
      )
    ) %>%
    make.epoch.labels(epoch.config) %>%
    select(
      generation,
      generation.start.year,
      generation.end.year,
      generation.window.years,
      aggregation.start.year,
      aggregation.end.year,
      hacker.window.years,
      is.truncated,
      population.boundary.year,
      periods,
      epoch,
      imported,
      births,
      birth.rate,
      pop.end,
      previous.birth.periods,
      births.in.previous.25.years,
      population.at.generation.start,
      female.total,
      male.total,
      adult.population.age25.proxy,
      parental.population,
      parental.population.method,
      k,
      ne
    )

  validate.analysis.outputs(
    AA.ne = AA.ne,
    AA.overlap = AA.overlap,
    AA.mixing.props = AA.mixing.props,
    epoch.config = epoch.config,
    hacker.start.year = hacker.start.year,
    gen.time = gen.time,
    present.year = present.year
  )

  list(
    AA.ne = AA.ne,
    AA.overlap = AA.overlap,
    AA.mixing.props = AA.mixing.props,
    epoch.config = epoch.config,
    ne.g1 = ne.g1,
    epoch.3.anchor.ne = epoch.3.anchor.ne,
    epoch.3.growth.rate = epoch.3.growth.rate
  )
}


# tunable parameters ----
gen.time <- 25
hacker.start.year <- 1610
hacker.end.year <- 1865
present.year <- 2000
epoch.1.len <- 8
epoch.2.len <- 3
epoch.3.len <- 5

hacker.tbl <- data.frame(
  period = c(
    "1610-1620", "1620-1630", "1630-1640", "1640-1650", "1650-1660",
    "1660-1670", "1670-1680", "1680-1690", "1690-1700", "1700-1710",
    "1710-1720", "1720-1730", "1730-1740", "1740-1750", "1750-1760",
    "1760-1770", "1770-1780", "1780-1790", "1790-1800", "1800-1810",
    "1810-1820", "1820-1830", "1830-1840", "1840-1850", "1850-1860",
    "1860-1861", "1861-1862", "1862-1863", "1863-1864", "1864-1865"
  ),
  pop.beginning = c(
    0, 20, 59, 585, 1568,
    2862, 4444, 6832, 16394, 27806,
    41844, 67294, 95669, 156040, 247027,
    319290, 459446, 558921, 706514, 908036,
    1195182, 1550757, 2021968, 2530405, 3204420,
    3953760, 3991133, 3903859, 3815760, 3726828
  ),
  pop.end = c(
    20, 59, 585, 1568, 2862,
    4444, 6832, 16394, 27806, 41844,
    67294, 95669, 156040, 247027, 319290,
    459446, 558921, 706514, 908036, 1195182,
    1550757, 2021968, 2530405, 3204420, 3953760,
    3991133, 3903859, 3815760, 3726828, 3630336
  ),
  growth.rate.10yr = c(
    NA, 195.0, 891.5, 168.0, 82.5,
    55.3, 53.7, 140.0, 69.6, 50.5,
    60.8, 42.2, 63.1, 58.3, 29.3,
    43.9, 21.7, 26.4, 28.5, 31.6,
    29.8, 30.4, 25.1, 26.6, 23.4,
    NA, NA, NA, NA, NA
  ),
  person.years.lived = c(
    16, 361, 2293, 9970, 21505,
    35952, 55527, 109243, 215999, 343482,
    535651, 806513, 1234030, 1980625, 2816150,
    3851269, 5075599, 6298380, 8030652, 10450423,
    13652609, 17759559, 22666906, 28541607, 35659777,
    3972417, 3947335, 3859642, 3771119, 3144251
  ),
  hours.worked.millions = c(
    0, 1, 5, 23, 49,
    82, 127, 251, 496, 788,
    1229, 1850, 2831, 4544, 6461,
    8835, 11644, 14449, 18423, 23975,
    31321, 40743, 52001, 65478, 81808,
    9113, 9056, 8855, 8651, 7213
  ),
  assumed.birth.rate = c(
    NA, 22.4, 22.4, 22.4, 22.4,
    22.4, 26.5, 27.3, 29.1, 31.0,
    34.6, 41.3, 41.4, 41.2, 46.4,
    39.5, 39.1, 44.7, 50.3, 55.8,
    55.7, 56.5, 56.5, 54.0, 52.3,
    49.6, 49.6, 49.6, 49.6, 49.6
  ),
  births = c(
    NA, 8, 51, 223, 481,
    804, 1472, 2986, 6295, 10635,
    18524, 33316, 51150, 81536, 130654,
    151957, 198660, 281549, 403648, 583395,
    759768, 1003859, 1280680, 1539820, 1865898,
    197032, 195788, 191438, 187048, 155955
  ),
  imported = c(
    32, 150, 110, 301, 1335,
    3038, 3582, 4947, 8197, 16442,
    17622, 41334, 66672, 26455, 44718,
    59557, 41129, 25105, 16183, 77368,
    5207, 410, 91, 0, 303,
    110, 0, 0, 0, 0
  ),
  smuggled = c(
    NA, NA, NA, NA, NA,
    NA, NA, NA, NA, NA,
    NA, NA, NA, NA, NA,
    NA, NA, NA, NA, NA,
    NA, 1100, 900, 500, 200,
    NA, NA, NA, NA, NA
  ),
  natural.increase.rate = c(
    NA, -30.8, 18.1, 6.8, -0.2,
    -4.0, -2.2, 4.2, 1.5, -0.7,
    1.5, -1.6, -0.5, 3.3, 1.0,
    2.1, 1.3, 2.2, 2.6, 2.4,
    2.7, 2.9, 2.4, 2.4, 2.2,
    NA, NA, NA, NA, NA
  ),
  births.plus.imports = c(
    32, 158, 161, 524, 1816,
    3842, 5054, 7933, 14492, 27077,
    36146, 74650, 117822, 107991, 175372,
    211514, 239789, 306654, 419831, 660763,
    764975, 1005369, 1281671, 1540320, 1866401,
    197142, 195788, 191438, 187048, 155955
  ),
  percent.total = c(
    0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.1, 0.1, 0.1, 0.3,
    0.4, 0.8, 1.2, 1.1, 1.8,
    2.2, 2.4, 3.1, 4.3, 6.7,
    7.8, 10.3, 13.1, 15.7, 19.0,
    2.0, 2.0, 2.0, 1.9, 1.6
  ),
  cumulative.number = c(
    32, 190, 351, 875, 2691,
    6533, 11587, 19520, 34012, 61089,
    97235, 171885, 289707, 397698, 573070,
    784584, 1024373, 1331027, 1750858, 2411621,
    3176596, 4181965, 5463636, 7003956, 8870357,
    9067499, 9263287, 9454725, 9641773, 9797728
  ),
  cumulative.percent = c(
    0.0, 0.0, 0.0, 0.0, 0.0,
    0.1, 0.1, 0.2, 0.3, 0.6,
    1.0, 1.8, 3.0, 4.1, 5.8,
    8.0, 10.5, 13.6, 17.9, 24.6,
    32.4, 42.7, 55.8, 71.5, 90.5,
    92.5, 94.5, 96.5, 98.4, 100.0
  ),
  stringsAsFactors = FALSE
)

afr.founder.count <- 183
afr.female.proportion <- 1 / 3
eur.founder.count <- 25
eur.female.proportion <- 0.5
epoch.3.growth.rate <- 0.023175

mooney.cg <- c(1, 1, 0.9835, 0.8602, 0.8551, 0.7826, 0.5380, 0.1418)
epoch.1.eur.prop <- 0.08
epoch.2.afr.prop <- 0.06
epoch.2.eur.prop <- 0.03
epoch.3.afr.prop <- 0
epoch.3.eur.prop <- 0


# calculations ----
analysis.results <- calc.aa.ne.results()

AA.ne <- analysis.results$AA.ne
AA.overlap <- analysis.results$AA.overlap
AA.mixing.props <- analysis.results$AA.mixing.props
epoch.config <- analysis.results$epoch.config
ne.g1 <- analysis.results$ne.g1
epoch.3.anchor.ne <- analysis.results$epoch.3.anchor.ne


# plotting ----
summary.generations <- c(
  epoch.config$epoch.1.generations,
  epoch.config$epoch.2.generations,
  epoch.config$final.generation
)
summary.generations <- sort(unique(summary.generations))

cat("\nHistorical generation-level ADMIX Ne summary\n")
print(
  AA.ne %>%
    filter(generation %in% summary.generations) %>%
    select(
      generation,
      generation.start.year,
      generation.end.year,
      epoch,
      is.truncated,
      imported,
      births,
      birth.rate,
      female.total,
      male.total,
      adult.population.age25.proxy,
      parental.population,
      parental.population.method,
      k,
      ne
    ),
  n = Inf
)

cat(
  paste0(
    "\nHistorical present-day G",
    epoch.config$final.generation,
    " Ne: ",
    round(AA.ne$ne[AA.ne$generation == epoch.config$final.generation], 1),
    "\n"
  )
)

AA.ne.plot <- ggplot(
  data = AA.ne,
  aes(x = generation.end.year, y = ne)
) +
  geom_line(color = "#386cb0", linewidth = 1.1) +
  geom_point(aes(shape = epoch), size = 2.5) +
  scale_y_continuous(labels = scales::comma) +
  scale_x_continuous(
    breaks = AA.ne$generation.end.year,
    labels = AA.ne$generation
  ) +
  labs(
    x = "Generation",
    y = "Effective population size",
    shape = NULL,
    title = "Historical generation-level African American Ne",
    subtitle = paste0(
      "G1 founder Ne = ",
      round(ne.g1, 1),
      " | G",
      epoch.config$epoch.3.anchor.generation,
      " anchor Ne = ",
      round(epoch.3.anchor.ne, 1),
      " | G",
      epoch.config$final.generation,
      " Ne = ",
      round(AA.ne$ne[AA.ne$generation == epoch.config$final.generation], 1)
    )
  ) +
  theme_bw()

if (interactive()) {
  print(AA.ne.plot)
}


# tennessen afr/eur ne trajectory ----
# constants from `other_scripts/const.sh` and `build_3T_demography.py`
t.af.years <- 148000
t.ooa.years <- 51000
t.eu0.years <- 23000
t.eg.years <- 5115
r.eu0 <- 0.00307
r.eu <- 0.0195
r.afr <- 0.0166
n.a <- 7310
n.afr.1 <- 14474
n.b <- 1861
n.eu0 <- 1032

t.af <- t.af.years / gen.time
t.ooa <- t.ooa.years / gen.time
t.eu0 <- t.eu0.years / gen.time
t.eg <- t.eg.years / gen.time
n.eu1 <- n.eu0 / exp(-r.eu0 * (t.eu0 - t.eg))
n.eu.present <- n.eu1 / exp(-r.eu * t.eg)
n.afr.present <- n.afr.1 / exp(-r.afr * t.eg)

tennessen.start.year <- min(AA.ne$generation.start.year)
tennessen.end.year <- present.year
tennessen.years <- seq(tennessen.start.year, tennessen.end.year, by = 1)
tennessen.time.bp <- tennessen.end.year - tennessen.years
tennessen.time.gen <- tennessen.time.bp / gen.time

tennessen.ne <- tibble(
  year = tennessen.years,
  time.bp = tennessen.time.bp,
  time.gen = tennessen.time.gen
) %>%
  mutate(
    afr = case_when(
      time.gen <= t.eg ~ n.afr.present * exp(-r.afr * time.gen),
      time.gen <= t.ooa ~ n.afr.1,
      time.gen <= t.af ~ n.a,
      TRUE ~ NA_real_
    ),
    eur = case_when(
      time.gen <= t.eg ~ n.eu.present * exp(-r.eu * time.gen),
      time.gen <= t.eu0 ~ n.eu1 * exp(-r.eu0 * (time.gen - t.eg)),
      time.gen <= t.ooa ~ n.b,
      TRUE ~ NA_real_
    )
  )

tennessen.ne.long <- tennessen.ne %>%
  pivot_longer(
    cols = c(afr, eur),
    names_to = "population",
    values_to = "ne"
  )

tennessen.ne.plot <- ggplot(
  tennessen.ne.long,
  aes(x = year, y = ne, color = population)
) +
  geom_line(linewidth = 1.1) +
  scale_y_continuous(labels = scales::comma) +
  scale_x_continuous(
    breaks = seq(tennessen.start.year, tennessen.end.year, by = 50)
  ) +
  labs(
    x = "Year",
    y = "Effective population size",
    color = NULL,
    title = paste0(
      "Tennessen AFR/EUR effective population size, ",
      tennessen.start.year,
      "-",
      tennessen.end.year
    )
  ) +
  theme_bw()

if (interactive()) {
  print(tennessen.ne.plot)
}
