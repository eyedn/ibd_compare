# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology 
# Mooney Lab
# ___
# kinship.mat2tbl.R
# ______________________________________________________________________________


library(dplyr)

kinship.mat2tbl <- function(mat, rep.idx, pop.name){
  
  # grab upper triangle only (no diagonal)
  idx <- which(upper.tri(mat), arr.ind = TRUE)

  if (!is.null(rownames(mat)) && !is.null(colnames(mat))) {
    sam1 <- rownames(mat)[idx[, 1]]
    sam2 <- colnames(mat)[idx[, 2]]
  } else {
    sam1 <- idx[, 1] - 1
    sam2 <- idx[, 2] - 1
  }
  
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
