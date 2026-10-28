# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# read.rep.calls.R
# ______________________________________________________________________________


source(if (file.exists("analysis/R_utils/read.figure.inputs.R")) {
  "analysis/R_utils/read.figure.inputs.R"
} else "R_utils/read.figure.inputs.R", local = TRUE)

read.rep.calls <- function(rep.idx, dir.paths.map, unrels.tbl, cores) {
  read.figure.calls(dir.paths.map, unrels.tbl, cores = cores,
                    rep.idx = rep.idx)
}
