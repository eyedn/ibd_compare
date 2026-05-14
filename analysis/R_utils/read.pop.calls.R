# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology 
# Mooney Lab
# ___
# read.pop.calls.R
# ______________________________________________________________________________


library(tidyverse)
library(furrr)
library(future)
library(progressr)

read.pop.calls <- function(focal.pop, dir.paths.map, unrels.tbl, cores) {
  # define table of files specific to one pop
  all.dir.vec <- names(dir.paths.map)
  file.tbl <- tibble(dir = all.dir.vec) %>%
    mutate(
      file = map(dir, ~ list.files(
        .x, full.names = TRUE,
        pattern = "_\\d+(\\.(hot|ibdends|hotibdends|ibdendshot))?\\.bed\\.gz$"
      ))
    ) %>%
    unnest(file) %>%
    filter(!is.na(file))
  file.tbl <- file.tbl %>%
    filter(str_detect(file, focal.pop))
  
  # parallel read + annotate per file
  plan(multisession, workers = cores)
  all.calls <- NA
  all.calls <- with_progress({
    p <- progressor(along = seq_len(nrow(file.tbl)))
    
    future_map_dfr(
      seq_len(nrow(file.tbl)),
      function(i) {
        p()
        
        dir <- file.tbl$dir[i]
        f   <- file.tbl$file[i]
        fname <- basename(f)
        parts <- strsplit(fname, "_")[[1]]
        
        if (parts[1] == "true") {
          imap <- "true"
          software <- "true"
        } else {
          imap <- parts[3]
          software <- parts[1]
        }
        
        pop <- parts[length(parts) - 1]
        rep <- as.integer(sub("(\\d+).*", "\\1", parts[length(parts)]))
        filt <- sub(
          ".*_\\d+(?:\\.(hot|ibdends|hotibdends|ibdendshot))?\\.bed\\.gz$",
          "\\1",
          fname
        )
        if (filt == "") filt <- "none"
        tmp <- readr::read_tsv(
          f,
          col_names = TRUE,
          show_col_types = FALSE,
          progress = FALSE
        ) %>%
          mutate(
            sam1 = as.character(sam1),
            sam2 = as.character(sam2)
          )
        
        tmp$pop <- pop
        tmp$imap <- imap
        tmp$filt <- filt
        tmp$rep <- rep
        tmp$software <- software
        tmp$software.code <- unname(dir.paths.map[dir])
        sc <- tmp$software.code[1]
        tmp$param.code <- if (sc == "true") "true" else paste0(substr(sc, 1, 1), stringr::str_extract(sc, "\\d+$"))
        ids <- unrels.tbl$ids[unrels.tbl$rep == rep & unrels.tbl$pop == pop][[1]]
        ids <- as.character(ids)
        tmp$is.unrel <- if (is.null(ids)) FALSE else (tmp$sam1 %in% ids & tmp$sam2 %in% ids)
        
        tmp
      },
      .options = furrr_options(seed = TRUE)
    )
  })
  plan(sequential)
  return(all.calls)
}
