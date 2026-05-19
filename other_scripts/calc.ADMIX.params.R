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

# rate 1 ----
# data from Hacker 2021
AA.pop.1.all <- data.frame(
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

AA.pop.1 <- AA.pop.1.all %>%
  mutate(
    start.year = as.numeric(sub("-.*", "", period)),
    imported.prev = dplyr::lag(imported),
    pop.beginning.no.prev.imports = pop.end - imported.prev,
    std.year = start.year - min(start.year) - 10
  ) %>%
  filter(
    !is.na(pop.beginning.no.prev.imports),
    start.year <= 1850
    ) %>%
  select(start.year, std.year, pop.beginning.no.prev.imports) 

# fit linear model on log-transformed size and extract coefficients
log.lm.model.1 <- lm(
  log(AA.pop.1$pop.beginning.no.prev.imports) ~ AA.pop.1$std.year
  )
summary(log.lm.model.1)
r0.1 <- coef(log.lm.model.1)[2]
i0.1 <- exp(coef(log.lm.model.1)[1])

log.nls.model.1 <- nls(
  data = AA.pop.1, pop.beginning.no.prev.imports ~ a * exp(b * std.year), 
  start = list(a = i0.1, b = r0.1)
)
summary(log.nls.model.1)
i.1 <- 1.207e+04 # from summary
r.1 <- 2.531e-02 # from summary


# rate 2 ----
# from...
#   1. Time Almanac 2005,
#   2. Statistical Abstract of the United States: 2003
#   3. We, the American Blacks
AA.pop.2 <- data.frame(
  year = c(
    1850, 1860, 1900, 1910, 1920, 1930, 1940, 1950, 1960, 1970, 1980, 1990, 2000
    ), 
  size = c(
    3.6, 4.4, 8.8, 9.8, 10.5, 11.9, 12.9, 15.0, 18.9, 22.6, 26.5, 30.0, 34.6
    )
  )
AA.pop.2 <- AA.pop.2 %>% 
  mutate(
    std.year = year - min(year),
    size = size * 1e6
    )

# fit linear model on log-transformed size and extract coefficients
log.lm.model.2 <- lm(
  log(AA.pop.2$size) ~ AA.pop.2$std.year
)
summary(log.lm.model.2)
r0.2 <- coef(log.lm.model.2)[2]
i0.2 <- exp(coef(log.lm.model.2)[1])


log.nls.model.2 <- nls(
  data = AA.pop.2, size ~ a * exp(b * std.year), 
  start = list(a = i0.2, b = r0.2)
)
summary(log.nls.model.2)
i.2 <- 3.630e+06 # from summary
r.2 <- 1.506e-02 # from summary


# timings ----
# from Hacker, Baharian, and Mooney
g.switch <- 9 # at 1850, 
g.tot <- 15 # 15 generations from 2000 is 1925. approx start time of admix pop

# generation-scale timings ----
gen.time <- 25
admix.start.year <- 1625
epoch.switch.year <- 1850
present.year <- 2000
g.switch.gen <- (epoch.switch.year - admix.start.year) / gen.time
g.tot.gen <- (present.year - admix.start.year) / gen.time

# pop size estimates ----
N.0 <- AA.pop.1$pop.beginning.no.prev.imports[2]
N.1 <- N.0 * exp(g.switch * r.1)
N.2 <- N.1 * exp((g.tot - g.switch) * r.2)

# build observed data for plotting ----
AA.pop.plot.obs <- bind_rows(
  AA.pop.1 %>%
    transmute(
      year = start.year,
      size = pop.beginning.no.prev.imports,
      epoch = "Epoch 1: 1620-1850"
    ),
  AA.pop.2 %>%
    transmute(
      year = year,
      size = size,
      epoch = "Epoch 2: 1850-2000"
    )
)

# build fitted curves for plotting ----
AA.pop.plot.fit.1 <- tibble(
  year = seq(min(AA.pop.1$start.year), max(AA.pop.1$start.year), by = 1)
) %>%
  mutate(
    std.year = year - 1620,
    size = i.1 * exp(r.1 * std.year),
    epoch = "Epoch 1 fit"
  )

AA.pop.plot.fit.2 <- tibble(
  year = seq(min(AA.pop.2$year), max(AA.pop.2$year), by = 1)
) %>%
  mutate(
    std.year = year - 1850,
    size = i.2 * exp(r.2 * std.year),
    epoch = "Epoch 2 fit"
  )

AA.pop.plot.fit <- bind_rows(AA.pop.plot.fit.1, AA.pop.plot.fit.2)

# plot ----
ggplot() +
  geom_point(
    data = AA.pop.plot.obs,
    aes(x = year, y = size, color = epoch),
    size = 2,
    alpha = 0.8
  ) +
  geom_line(
    data = AA.pop.plot.fit,
    aes(x = year, y = size, color = epoch),
    linewidth = 1.1
  ) +
  scale_y_continuous(labels = scales::comma) +
  scale_color_manual(
    values = c(
      "Epoch 1: 1625-1850" = "#1b9e77",
      "Epoch 2: 1850-2000" = "#d95f02",
      "Epoch 1 fit" = "#1b9e77",
      "Epoch 2 fit" = "#d95f02"
    )
  ) +
  labs(
    x = "Year",
    y = "Population size",
    color = NULL,
    title = "African American population size across two growth epochs",
    subtitle = paste0(
      "Epoch 1: r = ", round(r.1, 5),
      " | Epoch 2: r = ", round(r.2, 5)
    )
  ) +
  theme_bw()

# plausable Ne curves ----
# build fitted curves for plotting ----
AA.pop.plot.fitNe.1 <- tibble(
  year = seq(min(AA.pop.1$start.year), max(AA.pop.1$start.year), by = 1)
) %>%
  mutate(
    std.year = year - 1620,
    size = 435 * exp(r.1 * std.year),
    epoch = "Epoch 1 fit"
  )

AA.pop.plot.fitNe.2 <- tibble(
  year = seq(min(AA.pop.2$year), max(AA.pop.2$year), by = 1)
) %>%
  mutate(
    std.year = year - 1850,
    size = as.numeric(AA.pop.plot.fitNe.1[nrow(AA.pop.plot.fitNe.1), "size"]) 
    * exp(r.2 * std.year),
    epoch = "Epoch 2 fit"
  )

AA.pop.plot.fitNe <- bind_rows(AA.pop.plot.fitNe.1, AA.pop.plot.fitNe.2)

# plot Ne ----
ggplot() +
  geom_line(
    data = AA.pop.plot.fitNe,
    aes(x = year, y = size, color = epoch),
    linewidth = 1.1
  ) +
  scale_y_continuous(labels = scales::comma) +
  scale_color_manual(
    values = c(
      "Epoch 1: 1625-1850" = "#1b9e77",
      "Epoch 2: 1850-2000" = "#d95f02",
      "Epoch 1 fit" = "#1b9e77",
      "Epoch 2 fit" = "#d95f02"
    )
  ) +
  labs(
    x = "Year",
    y = "Population size",
    color = NULL,
    title = "African American population size across two growth epochs",
    subtitle = paste0(
      "Epoch 1: Ne(1620) = ", 
      round(as.numeric(AA.pop.plot.fitNe[1, "size"]), 1), 
      ", r = ", round(r.1, 5),
      " | Epoch 2: Ne(1850) = ", 
      round(as.numeric(AA.pop.plot.fitNe[nrow(AA.pop.plot.fitNe.1), "size"]), 1),
      ", r = ", round(r.2, 5),
      ", Ne(2000) = ",
      round(as.numeric(AA.pop.plot.fitNe[nrow(AA.pop.plot.fitNe), "size"]), 1)
    )
  ) +
  theme_bw()

# generation-scale fits ----
AA.pop.1.gen <- AA.pop.1.all %>%
  mutate(
    start.year = as.numeric(sub("-.*", "", period)),
    imported.prev = dplyr::lag(imported),
    pop.beginning.no.prev.imports = pop.end - imported.prev,
    std.gen = (start.year - admix.start.year) / gen.time
  ) %>%
  filter(
    !is.na(pop.beginning.no.prev.imports),
    start.year >= admix.start.year,
    start.year <= epoch.switch.year
  ) %>%
  select(start.year, std.gen, pop.beginning.no.prev.imports)

AA.pop.2.gen <- AA.pop.2 %>%
  transmute(
    year = year,
    std.gen = (year - epoch.switch.year) / gen.time,
    size = size
  )

# generation-scale rate 1 ----
log.lm.model.1.gen <- lm(
  log(AA.pop.1.gen$pop.beginning.no.prev.imports) ~ AA.pop.1.gen$std.gen
)
summary(log.lm.model.1.gen)
r0.1.gen <- coef(log.lm.model.1.gen)[2]
i0.1.gen <- exp(coef(log.lm.model.1.gen)[1])

log.nls.model.1.gen <- nls(
  data = AA.pop.1.gen,
  pop.beginning.no.prev.imports ~ a * exp(b * std.gen),
  start = list(a = i0.1.gen, b = r0.1.gen)
)
summary(log.nls.model.1.gen)
i.1.gen <- coef(log.nls.model.1.gen)[["a"]]
r.1.gen <- coef(log.nls.model.1.gen)[["b"]]

# generation-scale rate 2 ----
log.lm.model.2.gen <- lm(log(AA.pop.2.gen$size) ~ AA.pop.2.gen$std.gen)
summary(log.lm.model.2.gen)
r0.2.gen <- coef(log.lm.model.2.gen)[2]
i0.2.gen <- exp(coef(log.lm.model.2.gen)[1])

log.nls.model.2.gen <- nls(
  data = AA.pop.2.gen,
  size ~ a * exp(b * std.gen),
  start = list(a = i0.2.gen, b = r0.2.gen)
)
summary(log.nls.model.2.gen)
i.2.gen <- coef(log.nls.model.2.gen)[["a"]]
r.2.gen <- coef(log.nls.model.2.gen)[["b"]]

# generation-scale pop size estimates ----
N.0.gen <- 435
N.1.gen <- N.0.gen * exp(g.switch.gen * r.1.gen)
N.2.gen <- N.1.gen * exp((g.tot.gen - g.switch.gen) * r.2.gen)

# console comparison output ----
cat(
  paste0(
    "\nYear vs generation growth-rate comparison\n",
    "Epoch 1 yearly r: ", round(r.1, 8),
    " | Epoch 1 generational r: ", round(r.1.gen, 8), "\n",
    "Epoch 2 yearly r: ", round(r.2, 8),
    " | Epoch 2 generational r: ", round(r.2.gen, 8), "\n",
    "g.switch.gen: ", g.switch.gen,
    " | g.tot.gen: ", g.tot.gen, "\n\n"
  )
)

# generation-scale census curves ----
AA.pop.plot.fit.1.gen <- tibble(
  year = seq(min(AA.pop.1.gen$start.year), max(AA.pop.1.gen$start.year), by = 1)
) %>%
  mutate(
    std.gen = (year - admix.start.year) / gen.time,
    size = i.1.gen * exp(r.1.gen * std.gen),
    epoch = "Epoch 1 fit gen"
  )

AA.pop.plot.fit.2.gen <- tibble(
  year = seq(min(AA.pop.2.gen$year), max(AA.pop.2.gen$year), by = 1)
) %>%
  mutate(
    std.gen = (year - epoch.switch.year) / gen.time,
    size = i.2.gen * exp(r.2.gen * std.gen),
    epoch = "Epoch 2 fit gen"
  )

AA.pop.plot.fit.gen <- bind_rows(AA.pop.plot.fit.1.gen, AA.pop.plot.fit.2.gen)

# plot generation-scale census fits ----
ggplot() +
  geom_point(
    data = AA.pop.plot.obs,
    aes(x = year, y = size, color = epoch),
    size = 2,
    alpha = 0.8
  ) +
  geom_line(
    data = AA.pop.plot.fit.gen,
    aes(x = year, y = size, color = epoch),
    linewidth = 1.1
  ) +
  scale_y_continuous(labels = scales::comma) +
  scale_color_manual(
    values = c(
      "Epoch 1: 1625-1850" = "#1b9e77",
      "Epoch 2: 1850-2000" = "#d95f02",
      "Epoch 1: 1620-1850" = "#1b9e77",
      "Epoch 1 fit gen" = "#1b9e77",
      "Epoch 2 fit gen" = "#d95f02"
    )
  ) +
  labs(
    x = "Year",
    y = "Population size",
    color = NULL,
    title = "African American population size across two growth epochs",
    subtitle = paste0(
      "Epoch 1: r/gen = ", round(r.1.gen, 5),
      " | Epoch 2: r/gen = ", round(r.2.gen, 5)
    )
  ) +
  theme_bw()

# generation-scale plausable Ne curves ----
AA.pop.plot.fitNe.1.gen <- tibble(
  year = seq(min(AA.pop.1.gen$start.year), max(AA.pop.1.gen$start.year), by = 1)
) %>%
  mutate(
    std.gen = (year - admix.start.year) / gen.time,
    size = N.0.gen * exp(r.1.gen * std.gen),
    epoch = "Epoch 1 fit gen"
  )

AA.pop.plot.fitNe.2.gen <- tibble(
  year = seq(min(AA.pop.2.gen$year), max(AA.pop.2.gen$year), by = 1)
) %>%
  mutate(
    std.gen = (year - epoch.switch.year) / gen.time,
    size = N.1.gen * exp(r.2.gen * std.gen),
    epoch = "Epoch 2 fit gen"
  )

AA.pop.plot.fitNe.gen <- bind_rows(
  AA.pop.plot.fitNe.1.gen,
  AA.pop.plot.fitNe.2.gen
)

# plot generation-scale Ne ----
ggplot() +
  geom_line(
    data = AA.pop.plot.fitNe.gen,
    aes(x = year, y = size, color = epoch),
    linewidth = 1.1
  ) +
  scale_y_continuous(labels = scales::comma) +
  scale_color_manual(
    values = c(
      "Epoch 1 fit gen" = "#1b9e77",
      "Epoch 2 fit gen" = "#d95f02"
    )
  ) +
  labs(
    x = "Year",
    y = "Population size",
    color = NULL,
    title = "Generation-scale African American Ne across two growth epochs",
    subtitle = paste0(
      "Ne(1625) = ", round(N.0.gen, 1),
      ", r.1/gen = ", round(r.1.gen, 5),
      " | Ne(1850) = ", round(N.1.gen, 1),
      ", r.2/gen = ", round(r.2.gen, 5),
      ", Ne(2000) = ", round(N.2.gen, 1)
    )
  ) +
  theme_bw()

# historical generation-level Ne calculator ----
# Editable scenario inputs for the first historical Ne pass.
hist.afr.founder.count <- 326
hist.afr.female.proportion <- 0.5
hist.eur.founder.count <- 109
hist.eur.female.proportion <- 0.5
hist.epoch.3.growth.rate <- r.2.gen

# Stored for later msprime work. These proportions are not used below.
AA.pop.hist.mixing.proportions <- tibble(
  generation = 1:16,
  afr.proportion = rep(0.75, 16),
  eur.proportion = rep(0.25, 16)
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
    1610, 1630, 1660, 1700, 1750,
    1800, 1810, 1820, 1830, 1840
  ),
  generation.end.year = c(
    1630, 1660, 1700, 1750, 1800,
    1810, 1820, 1830, 1840, 1850
  ),
  periods = list(
    c("1610-1620", "1620-1630"),
    c("1630-1640", "1640-1650", "1650-1660"),
    c("1660-1670", "1670-1680", "1680-1690", "1690-1700"),
    c("1700-1710", "1710-1720", "1720-1730", "1730-1740",
      "1740-1750"),
    c("1750-1760", "1760-1770", "1770-1780", "1780-1790",
      "1790-1800"),
    c("1800-1810"),
    c("1810-1820"),
    c("1820-1830"),
    c("1830-1840"),
    c("1840-1850")
  )
)

AA.pop.hist.hacker <- AA.pop.1.all %>%
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
