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
  if (!file.exists(path)) stop("Missing required KING input: ", path)
  data <- readr::read_tsv(path, show_col_types = FALSE,
                          col_types = cols(.default = col_character()))
  names(data) <- sub("^#", "", names(data))
  required <- c("IID1", "IID2", "KINSHIP")
  if (!all(required %in% names(data))) {
    stop("Missing required KING columns in ", path, ": ",
         paste(setdiff(required, names(data)), collapse = ", "))
  }
  data %>% transmute(sam1 = IID1, sam2 = IID2,
                     kinship = as.numeric(KINSHIP))
}
