# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# calc.ADMIX.params.R
# ______________________________________________________________________________


library(tidyverse)

# historical generation-level Ne calculator ----
gen.time <- 25
admix.start.year <- 1625
epoch.switch.year <- 1850
present.year <- 2000

hacker2020.tbl <- data.frame(
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

# Editable scenario inputs for the first historical Ne pass.
hist.afr.founder.count <- 183
hist.afr.female.proportion <- 1/3
hist.eur.founder.count <- 25
hist.eur.female.proportion <- 0.5
hist.epoch.3.growth.rate <- 0.023175

# Stored for later msprime work. These proportions are not used below.
AA.pop.hist.mixing.proportions <- tibble(
  g = 1:16,
  s.afr = c(
    0.92, 0.90482, 0.791384, 0.786692, 0.720084, 
    0.06, 0.06, 0.06, 0.06, 0.06, 
    0, 0, 0, 0, 0, 0
  ),
  s.eur = c(
    0.08, 0.08, 0.08, 0.08, 0.08, 
    0.03, 0.03, 0.03, 0.03, 0.03, 
    0, 0, 0, 0, 0, 0
  ),
  h = c(
    0, 0.01518, 0.128616, 0.133308, 0.199916, 
    0.91, 0.91, 0.91, 0.91, 0.91, 
    0, 0, 0, 0, 0, 0
  )
)

calc.founder.ne <- function(founder.count, female.proportion) {
  female.count <- founder.count * female.proportion
  male.count <- founder.count - female.count
  
  if (female.count <= 0 || male.count <= 0) {
    return(0)
  }
  
  4 * female.count * male.count / (female.count + male.count)
}

AA.pop.hist.generation.map <- tibble(
  generation = 1:10,
  generation.start.year = c(
    1600, 1625, 1650, 1675, 1700,
    1725, 1750, 1775, 1800, 1825
  ),
  generation.end.year = c(
    1625, 1650, 1675, 1700, 1725,
    1750, 1775, 1800, 1825, 1850
  ),
  periods = list(
    c("1610-1620", "1620-1630"),
    c("1630-1640", "1640-1650"),
    c( "1650-1660", "1660-1670", "1670-1680"),
    c("1680-1690", "1690-1700"),
    c("1700-1710", "1710-1720", "1720-1730"), 
    c("1730-1740", "1740-1750"),
    c("1750-1760", "1760-1770", "1770-1780"), 
    c("1780-1790", "1790-1800"),
    c("1800-1810", "1810-1820", "1820-1830"),
    c("1830-1840", "1840-1850")
  )
)

AA.pop.hist.hacker <- hacker2020.tbl %>%
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
    pop.end
  )

AA.pop.hist.agg <- AA.pop.hist.generation.map %>%
  unnest(periods) %>%
  left_join(AA.pop.hist.hacker, by = c("periods" = "period")) %>%
  arrange(generation, start.year) %>%
  group_by(generation, generation.start.year, generation.end.year) %>%
  summarise(
    periods = paste(periods, collapse = ", "),
    imported = sum(imported, na.rm = TRUE),
    births = sum(births, na.rm = TRUE),
    birth.rate = weighted.mean(
      assumed.birth.rate,
      person.years.lived,
      na.rm = TRUE
    ),
    pop.end = dplyr::last(pop.end),
    .groups = "drop"
  ) %>%
  mutate(
    epoch = case_when(
      generation == 1 ~ "G1 founders",
      generation <= 5 ~ "G2-G5 imports plus births",
      TRUE ~ "G6-G10 birth-rate Ne"
    )
  )

AA.pop.hist.future <- tibble(
  generation = 11:16,
  generation.start.year = epoch.switch.year + ((generation - 11) * gen.time),
  generation.end.year = epoch.switch.year + ((generation - 10) * gen.time),
  periods = NA_character_,
  imported = NA_real_,
  births = NA_real_,
  birth.rate = NA_real_,
  pop.end = NA_real_,
  epoch = "G11-G16 epoch-3 growth"
)

AA.pop.hist.ne <- bind_rows(AA.pop.hist.agg, AA.pop.hist.future) %>%
  arrange(generation) %>%
  mutate(
    prev.births = dplyr::lag(births),
    female.imported = if_else(
      generation >= 2 & generation <= 5,
      ceiling(imported / 3),
      NA_real_
    ),
    male.imported = if_else(
      generation >= 2 & generation <= 5,
      floor(2 * imported / 3),
      NA_real_
    ),
    female.potential = if_else(
      generation >= 2 & generation <= 5,
      0.5 * prev.births,
      NA_real_
    ),
    male.potential = if_else(
      generation >= 2 & generation <= 5,
      0.5 * prev.births,
      NA_real_
    ),
    female.total = female.imported + female.potential,
    male.total = male.imported + male.potential,
    adult.population = if_else(
      generation >= 6 & generation <= 10,
      pop.end - births,
      NA_real_
    ),
    k = if_else(
      generation >= 6 & generation <= 10,
      gen.time * birth.rate / 1000,
      NA_real_
    )
  )

hist.ne.g1 <- calc.founder.ne(
  hist.afr.founder.count,
  hist.afr.female.proportion
) +
  calc.founder.ne(
    hist.eur.founder.count,
    hist.eur.female.proportion
  )

AA.pop.hist.ne <- AA.pop.hist.ne %>%
  mutate(
    ne = case_when(
      generation == 1 ~ hist.ne.g1,
      generation >= 2 & generation <= 5 ~
        4 * female.total * male.total / (female.total + male.total),
      generation >= 6 & generation <= 10 ~
        (k * adult.population - 1) / k,
      TRUE ~ NA_real_
    )
  )

hist.ne.g10 <- AA.pop.hist.ne %>%
  filter(generation == 10) %>%
  pull(ne)

AA.pop.hist.ne <- AA.pop.hist.ne %>%
  mutate(
    ne = if_else(
      generation >= 11,
      hist.ne.g10 * exp(hist.epoch.3.growth.rate * (generation - 10)),
      ne
    )
  ) %>%
  select(
    generation,
    generation.start.year,
    generation.end.year,
    periods,
    epoch,
    imported,
    births,
    birth.rate,
    pop.end,
    female.total,
    male.total,
    adult.population,
    k,
    ne
  )

stopifnot(
  identical(sort(AA.pop.hist.ne$generation), 1:16),
  all(!is.na(AA.pop.hist.ne$female.total[2:5])),
  all(!is.na(AA.pop.hist.ne$male.total[2:5])),
  all(!is.na(AA.pop.hist.ne$adult.population[6:10])),
  all(!is.na(AA.pop.hist.ne$k[6:10])),
  all(is.finite(AA.pop.hist.ne$ne)),
  all(AA.pop.hist.ne$ne > 0)
)

cat("\nHistorical generation-level ADMIX Ne summary\n")
print(
  AA.pop.hist.ne %>%
    filter(generation <= 10 | generation == 16) %>%
    select(
      generation,
      generation.end.year,
      epoch,
      imported,
      births,
      birth.rate,
      female.total,
      male.total,
      adult.population,
      k,
      ne
    ),
  n = Inf
)

cat(
  paste0(
    "\nHistorical present-day G16 Ne: ",
    round(AA.pop.hist.ne$ne[AA.pop.hist.ne$generation == 16], 1),
    "\n"
  )
)

ggplot(
  data = AA.pop.hist.ne,
  aes(x = generation.end.year, y = ne)
) +
  geom_line(color = "#386cb0", linewidth = 1.1) +
  geom_point(aes(shape = epoch), size = 2.5) +
  scale_y_continuous(labels = scales::comma) +
  scale_x_continuous(
    breaks = AA.pop.hist.ne$generation.end.year,
    labels = AA.pop.hist.ne$generation
  ) +
  labs(
    x = "Generation",
    y = "Effective population size",
    shape = NULL,
    title = "Historical generation-level African American Ne",
    subtitle = paste0(
      "G1 founder Ne = ", round(hist.ne.g1, 1),
      " | G10 Ne = ", round(hist.ne.g10, 1),
      " | G16 Ne = ",
      round(AA.pop.hist.ne$ne[AA.pop.hist.ne$generation == 16], 1)
    )
  ) +
  theme_bw()
