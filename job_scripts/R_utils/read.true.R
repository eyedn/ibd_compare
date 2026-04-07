# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology 
# Mooney Lab
# ___
# read.true.R
# ______________________________________________________________________________


library(tidyverse)
source("job_scripts/R_utils/interp.cm.R")

read.true <- function(true.fname,  min.bp.len = 0, min.cM.len = 0, gmap) {
  # read in true segments
  true.segs <- readr::read_tsv(
    true.fname, col_names = FALSE, show_col_types = FALSE
    )
  true.segs <- true.segs %>%
    rename(
      chr = X1,
      bp.str = X2,
      bp.end = X3,
      sam1 = X4,
      hap1 = X5,
      sam2 = X6,
      hap2 = X7,
      node.mrca = X8,
      time.mrca = X9,
      bp.len = X10
    ) %>%
    mutate(
      bp.str = as.integer(bp.str),
      bp.end = as.integer(bp.end),
      sam1 = as.integer(str_extract(sam1, "\\d+$")),
      hap1   = as.integer(hap1),
      sam2 = as.integer(str_extract(sam2, "\\d+$")),,
      hap2   = as.integer(hap2),
      bp.len = as.integer(bp.len)
    ) %>% 
    select(-node.mrca, -time.mrca)
  print(paste0("total true segments: ", nrow(true.segs)))

  # read in gmap
  true.gmap <- readr::read_table(
    file = gmap, col_names = FALSE, show_col_types = FALSE
  )
  true.gmap <- true.gmap %>% rename(chr = X1, snp = X2, cM = X3, bp = X4)
  
  # generate cM.str and cM.end for each ibd call using linear interpolation
  true.segs <- true.segs %>%
    mutate(
      cM.str = interp.cm(bp.str, true.gmap$bp, true.gmap$cM),
      cM.end = interp.cm(bp.end, true.gmap$bp, true.gmap$cM),
      cM.len = cM.end - cM.str
    )
  
  # finalize standard format for future comparison
  true.segs <- true.segs %>%
    filter(bp.len >= min.bp.len, cM.len >= min.cM.len)
  print(paste0("total filtered segments: ", nrow(true.segs)))
  true.segs %>%
    # enforce sam1 <= sam2 lexicographically, swapping haps as needed
    mutate(
      swap   = sam1 > sam2,                          # TRUE if we need to swap
      sam1_n = if_else(swap, sam2, sam1),
      sam2_n = if_else(swap, sam1, sam2),
      hap1_n = if_else(swap, hap2, hap1),
      hap2_n = if_else(swap, hap1, hap2)
    ) %>%
    select(-sam1, -sam2, -hap1, -hap2, -swap) %>%    # drop old + helper
    rename(
      sam1 = sam1_n,
      sam2 = sam2_n,
      hap1 = hap1_n,
      hap2 = hap2_n
    ) %>%
    # order like a BED file
    arrange(bp.str, bp.end) %>%
    # select columns in desired order
    select(
      chr,
      bp.str,
      bp.end,
      bp.len,
      cM.str,
      cM.end,
      cM.len,
      sam1,
      hap1,
      sam2,
      hap2
    )
}
