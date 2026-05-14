# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology 
# Mooney Lab
# ___
# compute.kinship.coeffs.R
# ______________________________________________________________________________


library(tidyverse)

compute.kinship.coeffs <- function(ibd, n.diploids, genome.seg.len.cM,
                                   sample.id.mode = c(
                                     "simulation", "empirical"
                                   ),
                                   sample.ids = NULL) {
  sample.id.mode <- match.arg(sample.id.mode)
  
  message("computing IBD-based kinship with ", n.diploids, " diploids over a ", genome.seg.len.cM, "cM length of genome")
  req.cols <- c("chr", "sam1", "hap1", "sam2", "hap2",
                "bp.str", "bp.end", "cM.str", "cM.end")
  
  if (!all(req.cols %in% names(ibd))) {
    stop("ibd must contain columns: ", paste(req.cols, collapse = ", "))
  }

  if (sample.id.mode == "simulation") {
    if (is.null(n.diploids)) {
      stop("simulation kinship requires n.diploids")
    }
  } else {
    if (is.null(sample.ids)) {
      sample.ids <- sort(unique(c(ibd$sam1, ibd$sam2)))
    }
    sample.ids <- as.character(sample.ids)
    n.diploids <- length(sample.ids)
  }

  # extract necessary coli,ms
  ibd.condense <- ibd %>%
    select(chr, sam1, hap1, sam2, hap2, cM.str, cM.end) %>%
    mutate(
      sam1 = as.character(sam1),
      sam2 = as.character(sam2),
      sam.pair = paste0(sam1, "-", sam2)
    ) %>%
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
  if (sample.id.mode == "simulation") {
    kin.mat <- matrix(0, nrow = n.diploids, ncol = n.diploids)
  } else {
    kin.mat <- matrix(
      0,
      nrow = n.diploids,
      ncol = n.diploids,
      dimnames = list(sample.ids, sample.ids)
    )
  }

  for (i in seq_len(nrow(pair.summary))) {
    if (sample.id.mode == "simulation") {
      s1 <- suppressWarnings(as.integer(pair.summary$sam1[i])) + 1L
      s2 <- suppressWarnings(as.integer(pair.summary$sam2[i])) + 1L

      if (is.na(s1) || is.na(s2)) next
      if (s1 < 1L || s1 > n.diploids || s2 < 1L || s2 > n.diploids) next
    } else {
      s1 <- match(pair.summary$sam1[i], sample.ids)
      s2 <- match(pair.summary$sam2[i], sample.ids)

      if (is.na(s1) || is.na(s2)) next
    }

    kin.mat[s1, s2] <- pair.summary$kinship[i]
    kin.mat[s2, s1] <- pair.summary$kinship[i]
  }
  
  diag(kin.mat) <- 0.5

  if (sample.id.mode == "simulation") {
    pair.summary <- pair.summary %>%
      mutate(
        sam1 = as.integer(sam1),
        sam2 = as.integer(sam2)
      )
  }
  
  list(
    kinship.long   = pair.summary,
    kinship.matrix = kin.mat
  )
}
