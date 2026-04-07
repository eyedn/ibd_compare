# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology 
# Mooney Lab
# ___
# read.rep.ibd.kinship.R
# ______________________________________________________________________________


library(tidyverse)
library(furrr)
library(future)
library(progressr)
source("R_utils/kinship.mat2tbl.R")

read.rep.ibd.kinship <- function(rep.idx, dir.paths.map, cores,
                                 which = c("inf", "true")) {
  
  which <- match.arg(which)
  all.dir.vec <- names(dir.paths.map)
  
  rep.file.tbl <- tibble(dir = all.dir.vec) %>%
    mutate(
      kin.dir = file.path(dir, "ibd_kinship_mats"),
      file = map(kin.dir, ~ list.files(
        .x, full.names = TRUE,
        pattern = paste0("^", which, "_[^_]+_", rep.idx, "\\.[^.]+\\.rds$")
      ))
    ) %>%
    unnest(file) %>%
    filter(!is.na(file)) %>%
    filter(!str_detect(basename(file), "ibdendshot"))  
  
  plan(multisession, workers = cores)
  
  kin.res <- with_progress({
    p <- progressor(along = seq_len(nrow(rep.file.tbl)))
    
    future_map_dfr(
      seq_len(nrow(rep.file.tbl)),
      function(i) {
        p()
        
        dir   <- rep.file.tbl$dir[i]
        f     <- rep.file.tbl$file[i]
        fname <- basename(f)
        
        # parse: {true|inf}_{pop}_{REP}.{filt}.rds
        m <- str_match(fname, "^(true|inf)_([^_]+)_(\\d+)\\.([^.]+)\\.rds$")
        pop  <- m[, 3]
        rep  <- as.integer(m[, 4])
        filt <- m[, 5]
        
        mat <- readRDS(f)
        kin.tbl <- kinship.mat2tbl(mat, rep, pop)
        
        software.code <- unname(dir.paths.map[dir])
        if (is.na(software.code) || length(software.code) == 0) software.code <- type
        
        kin.tbl %>%
          mutate(
            filt = filt,
            software.code = software.code
          )
      },
      .options = furrr_options(seed = TRUE)
    )
  })
  
  plan(sequential)
  kin.res
}
