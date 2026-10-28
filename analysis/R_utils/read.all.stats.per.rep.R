# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# read.all.stats.per.rep.R
# ______________________________________________________________________________


source(if (file.exists("analysis/R_utils/read.all.stats.R")) {
  "analysis/R_utils/read.all.stats.R"
} else "R_utils/read.all.stats.R", local = TRUE)

read.all.stats.per.rep <- function(file, meta.names, param.map) {
  data <- read.all.stats(file, param.map)
  missing <- setdiff(meta.names, names(data))
  if (length(missing)) {
    stop("Missing columns in ", file, ": ", paste(missing, collapse = ", "))
  }
  data
}
