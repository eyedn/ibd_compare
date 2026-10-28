# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# read.all.stats.R
# ______________________________________________________________________________


library(tidyverse)

read.all.stats <- function(file, param.map) {
  if (!file.exists(file)) stop("Missing required input: ", file)
  data <- readr::read_tsv(file, show_col_types = FALSE)
  # Exclude the retired configuration without renaming its results.
  data <- data %>% filter(is.na(software.with.params) | software.with.params !=
    "hapibd_minseed2.0_gap2000_extend1.0_out2.0_mark100_mac2")
  unknown <- setdiff(unique(data$software.with.params), names(param.map))
  if (length(unknown)) {
    stop("Unknown parameterization in ", file, ": ",
         paste(unknown, collapse = ", "))
  }
  data %>% mutate(software.code = unname(param.map[software.with.params]))
}
