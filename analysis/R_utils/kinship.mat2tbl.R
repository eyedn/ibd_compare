# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology 
# Mooney Lab
# ___
# kinship.mat2tbl.R
# ______________________________________________________________________________


# pattern: Functional Core

library(dplyr)

kinship.mat2tbl <- function(mat, rep.idx, pop.name){
  
  # Require global sample IDs on simulation kinship matrices.
  if (is.null(rownames(mat)) || is.null(colnames(mat))) {
    stop(
      "kinship matrices must have row and column names with global ",
      "sample IDs. Regenerate stage-4 ibd_kinship_mats outputs."
    )
  }

  # grab upper triangle only (no diagonal)
  idx <- which(upper.tri(mat), arr.ind = TRUE)

  sam1 <- rownames(mat)[idx[, 1]]
  sam2 <- colnames(mat)[idx[, 2]]

  tbl <- tibble(
    rep = rep.idx,
    pop = pop.name,
    sam1 = sam1,
    sam2 = sam2,
    kinship = mat[idx]
  )
  
  # create sam.pair column
  tbl %>%
    mutate(
      sam1.min = pmin(sam1, sam2),
      sam2.max = pmax(sam1, sam2),
      sam.pair = paste0(sam1.min, "-", sam2.max)
    ) %>%
    select(-sam1.min, -sam2.max)
}
