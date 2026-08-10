# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology 
# Mooney Lab
# ___
# read.plink2.king.R
# ______________________________________________________________________________


library(tidyverse)

read.plink2.king <- function(path) {
  readr::read_tsv(path, col_names = TRUE, show_col_types = FALSE) %>%
    rename(
      sam1    = `#IID1`,
      sam2    = IID2,
      kinship = KINSHIP
    ) %>%
    mutate(
      sam1 = as.character(sam1),
      sam2 = as.character(sam2)
    ) %>%
    select(sam1, sam2, kinship)
}
