# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# read.pop.calls.R
# ______________________________________________________________________________


source(if (file.exists("analysis/R_utils/read.figure.inputs.R")) {
  "analysis/R_utils/read.figure.inputs.R"
} else "R_utils/read.figure.inputs.R", local = TRUE)

read.pop.calls <- function(focal.pop, dir.paths.map, unrels.tbl, cores) {
  read.figure.calls(dir.paths.map, unrels.tbl, cores = cores,
                    focal.pop = focal.pop)
}
