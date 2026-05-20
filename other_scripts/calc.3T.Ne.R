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

# historical generation-level Ne calculator ----
gen.time <- 25
hacker.start.year <- 1610
hacker.end.year <- 1865
epoch.3.start.year <- 1875
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

# editable scenario inputs for the first historical Ne pass.
afr.founder.count <- 183
afr.female.proportion <- 1/3
eur.founder.count <- 25
eur.female.proportion <- 0.5
epoch.3.growth.rate <- 0.023175

# store for later msprime work; these proportions are not used below.
mooney.cg <- c(1, 1, 0.9835, 0.8602, 0.8551, 0.7826, 0.5380, 0.1418)
epoch.1.eur.prop <- 0.08
epoch.2.afr.prop <- 0.06
epoch.2.eur.prop <- 0.03
epoch.3.afr.prop <- 0
epoch.3.eur.prop <- 0
AA.mixing.props <- tibble(g = 1:16, s.afr = NA, s.eur = NA, h = NA)
AA.mixing.props <- AA.mixing.props %>%
  mutate(
    s.eur = c(
      rep(epoch.1.eur.prop, epoch.1.len), rep(epoch.2.eur.prop, epoch.2.len), 
      rep(epoch.3.eur.prop, epoch.3.len)
      ),
    s.afr = c(
      (1 - epoch.1.eur.prop) * mooney.cg, rep(epoch.2.afr.prop, epoch.2.len),
      rep(epoch.3.afr.prop, epoch.3.len)
      ),
    h = pmax(1 - s.afr - s.eur, 0)
    )

calc.founder.ne <- function(founder.count, female.proportion) {
  female.count <- founder.count * female.proportion
  male.count <- founder.count - female.count
  
  if (female.count <= 0 || male.count <= 0) {
    return(0)
  }
  
  4 * female.count * male.count / (female.count + male.count)
}

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

AA.generation.map <- tibble(
  generation = 1:16,
  generation.start.year = seq(1600, present.year - gen.time, by = gen.time),
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
      generation <= 11,
      pmin(generation.end.year, hacker.end.year),
      NA_real_
    ),
    pop.end = map_dbl(
      population.boundary.year,
      calc.boundary.population,
      hacker.tbl = AA.hacker.tbl
    ),
    is.truncated = !is.na(hacker.window.years) &
      hacker.window.years < generation.window.years,
    epoch = case_when(
      generation <= 10 ~
        "Epoch 1: 1600-1865 Hacker-table historical reconstruction",
      generation == 11 ~ "Epoch 2: 1865-1875 bridge period",
      TRUE ~ "Epoch 3: 1875-2000 exponential growth"
    )
  )

AA.parental.windows <- AA.generation.map %>%
  transmute(
    generation,
    generation.start.year,
    generation.end.year,
    window.start.year = generation.start.year - gen.time,
    window.end.year = generation.start.year
  )

AA.parental.births <- calc.period.overlaps(
  AA.parental.windows,
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
    births
  )

AA.parental.agg <- AA.parental.births %>%
  group_by(generation) %>%
  summarise(
    previous.birth.periods = paste(period, collapse = ", "),
    births.in.previous.25.years = sum(births * overlap.weight, na.rm = TRUE),
    .groups = "drop"
  )

AA.ne <- AA.agg %>%
  arrange(generation) %>%
  left_join(AA.parental.agg, by = "generation") %>%
  mutate(
    population.at.generation.start = map_dbl(
      generation.start.year,
      calc.boundary.population,
      hacker.tbl = AA.hacker.tbl
    ),
    prev.births = dplyr::lag(births),
    female.imported = if_else(
      generation >= 2 & generation <= epoch.1.len,
      ceiling(imported / 3),
      NA_real_
    ),
    male.imported = if_else(
      generation >= 2 & generation <= epoch.1.len,
      floor(2 * imported / 3),
      NA_real_
    ),
    female.potential = if_else(
      generation >= 2 & generation <= epoch.1.len,
      0.5 * prev.births,
      NA_real_
    ),
    male.potential = if_else(
      generation >= 2 & generation <= epoch.1.len,
      0.5 * prev.births,
      NA_real_
    ),
    female.total = female.imported + female.potential,
    male.total = male.imported + male.potential,
    adult.population.age25.proxy = if_else(
      generation > epoch.1.len & generation <= epoch.1.len + epoch.2.len,
      population.at.generation.start - births.in.previous.25.years,
      NA_real_
    ),
    parental.population = case_when(
      generation < epoch.1.len | generation > epoch.1.len + epoch.2.len ~ NA_real_,
      adult.population.age25.proxy > 0 ~ adult.population.age25.proxy,
      TRUE ~ population.at.generation.start
    ),
    parental.population.method = case_when(
      generation < epoch.1.len | generation > epoch.1.len + epoch.2.len ~ NA_character_,
      adult.population.age25.proxy > 0 ~ "age-25 proxy",
      TRUE ~ "generation-start population fallback"
    ),
    k = if_else(
      generation >= epoch.1.len & generation <= epoch.1.len + epoch.2.len,
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
      generation >= 2 & generation <= epoch.1.len ~
        4 * female.total * male.total / (female.total + male.total),
      generation > epoch.1.len & generation <= epoch.1.len + epoch.2.len ~
        (k * parental.population - 1) / k,
      TRUE ~ NA_real_
    )
  )

ne.g11 <- AA.ne %>%
  filter(generation == 11) %>%
  pull(ne)

AA.ne <- AA.ne %>%
  mutate(
    ne = if_else(
      generation.start.year >= epoch.3.start.year,
      ne.g11 * exp(epoch.3.growth.rate * (generation - 11)),
      ne
    )
  ) %>%
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

stopifnot(
  identical(sort(AA.ne$generation), 1:16),
  all(AA.ne$generation.window.years == gen.time),
  AA.ne$aggregation.start.year[1] == hacker.start.year,
  all(AA.overlap$overlap.weight > 0),
  all(AA.overlap$overlap.weight <= 1),
  all(!is.na(AA.ne$female.total[2:5])),
  all(!is.na(AA.ne$male.total[2:5])),
  all(!is.na(AA.ne$parental.population[6:11])),
  all(!is.na(AA.ne$k[6:11])),
  all(
    AA.ne$epoch[12:16] ==
      "Epoch 3: 1875-2000 exponential growth"
  ),
  all(is.finite(AA.ne$ne)),
  all(AA.ne$ne > 0)
)

cat("\nHistorical generation-level ADMIX Ne summary\n")
print(
  AA.ne %>%
    filter(generation <= epoch.1.len + epoch.2.len | generation == epoch.1.len + epoch.2.len + epoch.3.len) %>%
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
    "\nHistorical present-day G16 Ne: ",
    round(AA.ne$ne[AA.ne$generation == epoch.1.len + epoch.2.len + epoch.3.len], 1),
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
      "G1 founder Ne = ", round(ne.g1, 1),
      " | G11 Ne = ", round(ne.g11, 1),
      " | G16 Ne = ",
      round(AA.ne$ne[AA.ne$generation == 16], 1)
    )
  ) +
  theme_bw()

if (interactive()) {
  print(AA.ne.plot)
}
