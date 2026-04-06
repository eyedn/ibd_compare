# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology 
# Mooney Lab
# ___
# compute.kinship.coeffs.R ONE-BASED VER!!!
# ______________________________________________________________________________

library(tidyverse)

compute.kinship.coeffs <- function(ibd, n.diploids, genome.seg.len.cM) {
  
  message("computing IBD-based kinship with ", n.diploids, " diploids over a ", genome.seg.len.cM, "cM length of genome")
  req.cols <- c("chr", "sam1", "hap1", "sam2", "hap2",
                "bp.str", "bp.end", "cM.str", "cM.end")
  
  if (!all(req.cols %in% names(ibd))) {
    stop("ibd must contain columns: ", paste(req.cols, collapse = ", "))
  }

  # extract necessary coli,ms
  ibd.condense <- ibd %>%
    mutate(
      sam1 = as.integer(sam1),
      sam2 = as.integer(sam2)
    ) %>%
    select(chr, sam1, hap1, sam2, hap2, cM.str, cM.end) %>%
    mutate(sam.pair = paste0(as.character(sam1), "-", as.character(sam2))) %>%
    distinct()

  # summarize per diploid pair
  pair.summary <- ibd.condense %>%
    group_by(sam1, sam2) %>%
    group_modify(~ {
      df <- .x
      breaks <- sort(unique(c(df$cM.str, df$cM.end)))
      if (length(breaks) < 2) {
        return(tibble(ibd1 = 0, ibd2 = 0))
      }
      
      intervals <- tibble(
        cM.str = head(breaks, -1),
        cM.end = tail(breaks, -1),
        cM.len = cM.end - cM.str
      )
      
      counts <- intervals %>%
        rowwise() %>%
        mutate(
          n.sam.pairs = n_distinct(
            df$sam.pair[df$cM.str <= cM.str & df$cM.end >= cM.end]
            )
        ) %>%
        ungroup()
      
      tibble(
        ibd1 = sum(counts$cM.len[counts$n.sam.pairs == 1]),
        ibd2 = sum(counts$cM.len[counts$n.sam.pairs >= 2])
      )
    }) %>%
    ungroup()
  
  # kinship
  pair.summary <- pair.summary %>%
    mutate(
      kinship = (ibd1 + 2 * ibd2) / (4 * genome.seg.len.cM)
    )
  
  # matrix output
  kin.mat <- matrix(0, nrow = n.diploids, ncol = n.diploids)
  
  for (i in seq_len(nrow(pair.summary))) {
    s1 <- as.integer(pair.summary$sam1[i]) + 1L
    s2 <- as.integer(pair.summary$sam2[i]) + 1L
    
    # guard rails
    if (s1 < 1L || s1 > n.diploids || s2 < 1L || s2 > n.diploids) next
    
    kin.mat[s1, s2] <- pair.summary$kinship[i]
    kin.mat[s2, s1] <- pair.summary$kinship[i]
  }
  
  diag(kin.mat) <- 0.5
  
  list(
    kinship.long   = pair.summary,
    kinship.matrix = kin.mat
  )
}
