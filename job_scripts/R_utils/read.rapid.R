# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology 
# Mooney Lab
# ___
# read.rapid.R
# ______________________________________________________________________________


library(tidyverse)

read.rapid.bed <- function(fname, min.bp.len = 0, min.cM.len = 0, gmap) {
  # read in ibd calls
  rapid.calls <- readr::read_tsv(
    file = fname, col_names = FALSE, show_col_types = FALSE
  )
  print(paste0("total raw calls: ", nrow(rapid.calls)))
  
  rapid.calls <- rapid.calls %>%
    rename(chr = X1,
           bp.str = X2,
           bp.end = X3,
           sam1 = X4,
           hap1 = X5,
           sam2 = X6,
           hap2 = X7,
           cM.len = X8,
           idx.str = X9,
           idx.end = X10) %>%
    mutate(
      bp.str = as.integer(bp.str),
      bp.end = as.integer(bp.end),
      sam1 = as.integer(str_extract(sam1, "\\d+$")),
      hap1 = as.integer(hap1) + 1,
      sam2 = as.integer(str_extract(sam2, "\\d+$")),
      hap2 = as.integer(hap2) + 1,
      cM.len = as.numeric(cM.len),
      bp.len = bp.end - bp.str,
      idx.str = as.integer(idx.str),
      idx.end = as.integer(idx.end)
    ) %>%
    filter(bp.len >= min.bp.len, cM.len >= min.cM.len)
  print(paste0("total filtered calls: ", nrow(rapid.calls)))
  
  # read in gmap
  rapid.gmap <- readr::read_tsv(
    file = gmap, col_names = FALSE, show_col_types = FALSE
    )
  rapid.gmap <- rapid.gmap %>% rename(idx = X1, cM = X2)
  
  # generate cM.str and cM.end for each ibd call using the rapid map
  rapid.calls <- rapid.calls %>%
    left_join(rapid.gmap %>% rename(idx.str = idx, cM.str = cM), 
              by = "idx.str") %>%
    left_join(rapid.gmap %>% rename(idx.end = idx, cM.end = cM), 
              by = "idx.end")
  
  # finalize standard format for future comparison
  rapid.calls %>%
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
