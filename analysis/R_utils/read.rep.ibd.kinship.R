# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# read.rep.ibd.kinship.R
# ______________________________________________________________________________


source(if (file.exists("analysis/R_utils/read.figure.inputs.R")) {
  "analysis/R_utils/read.figure.inputs.R"
} else "R_utils/read.figure.inputs.R", local = TRUE)
source(if (file.exists("analysis/R_utils/kinship.mat2tbl.R")) {
  "analysis/R_utils/kinship.mat2tbl.R"
} else "R_utils/kinship.mat2tbl.R", local = TRUE)

read.rep.ibd.kinship <- function(rep.idx, dir.paths.map, cores,
    which = c("inf", "true"), true.dir = NULL) {
  read.figure.kinship(dir.paths.map, cores = cores, which = match.arg(which),
    true.dir = true.dir, rep.idx = rep.idx)
}
