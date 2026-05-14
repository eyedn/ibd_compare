# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology 
# Mooney Lab
# ___
# read.hapibd.R
# ______________________________________________________________________________


library(tidyverse)
source("job_scripts/R_utils/interp.cM.R")
source("job_scripts/R_utils/parse.sample.ids.R")

read.hapibd.bed <- function(fname, min.bp.len = 0, min.cM.len = 0, gmap,
                            sample.id.mode = c("simulation", "empirical")) {
  sample.id.mode <- match.arg(sample.id.mode)

  # read in ibd calls
  hapibd.calls <- readr::read_tsv(
    file = fname, 
    col_names = FALSE, show_col_types = FALSE)
  print(paste0("total raw calls: ", nrow(hapibd.calls)))
  
  hapibd.calls <- hapibd.calls %>%
    rename(chr = X1,
           bp.str = X2,
           bp.end = X3,
           sam1 = X4,
           hap1 = X5,
           sam2 = X6,
           hap2 = X7,
           cM.len = X8) %>%
    mutate(
      bp.str = as.integer(bp.str),
      bp.end = as.integer(bp.end),
      sam1 = parse.sample.ids(
        sam1,
        sample.id.mode,
        "hap-IBD simulation input contains non-integer sample IDs"
      ),
      hap1 = as.integer(hap1),
      sam2 = parse.sample.ids(
        sam2,
        sample.id.mode,
        "hap-IBD simulation input contains non-integer sample IDs"
      ),
      hap2 = as.integer(hap2),
      cM.len = as.numeric(cM.len),
      bp.len = bp.end - bp.str
    ) %>%
    filter(bp.len >= min.bp.len, cM.len >= min.cM.len)
  print(paste0("total filtered calls: ", nrow(hapibd.calls)))
  
  # read in gmap
  hapmap.gmap <- readr::read_table(
    file = gmap, col_names = FALSE, show_col_types = FALSE
  )
  hapmap.gmap <- hapmap.gmap %>% rename(chr = X1, snp = X2, cM = X3, bp = X4)
  
  # generate cM.str and cM.end for each ibd call using linear interpolation
  hapibd.calls <- hapibd.calls %>%
    mutate(
      cM.str = interp.cm(bp.str, hapmap.gmap$bp, hapmap.gmap$cM),
      cM.end = interp.cm(bp.end, hapmap.gmap$bp, hapmap.gmap$cM)
    )
  
  # finalize standard format for future comparison
  hapibd.calls %>%
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
