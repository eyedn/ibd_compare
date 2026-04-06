# ______________________________________________________________________________
# Aydin Karatas 
# ___
# University of Southern California 
# Department of Quantitative and Computational Biology 
# Mooney Lab 
# ___
# read.pop.ibd.kinship.R 
# ______________________________________________________________________________

library(tidyverse)
library(furrr)
library(future)
library(progressr)
source("R_utils/kinship.mat2tbl.R")

read.pop.ibd.kinship <- function(focal.pop, dir.paths.map, cores,
                                 which = c("true", "inf"),
                                 true.dir = NULL) {
  which <- match.arg(which)
  
  # choose dirs to search
  if (which == "true") {
    # truth is identical across methods/dirs -> only read once
    if (is.null(true.dir)) {
      # default: use the first directory in dir.paths.map
      true.dir <- names(dir.paths.map)[1]
    }
    search.dirs <- true.dir
    file.tbl <- tibble(dir = search.dirs) %>%
      mutate(
        kin.dir = file.path(dir, "ibd_kinship_mats"),
        file = map(kin.dir, ~ list.files(
          .x, full.names = TRUE,
          pattern = "\\.rds$"
        ))
      ) %>%
      unnest(file) %>%
      filter(!is.na(file)) %>%
      filter(str_detect(file, "none")) %>%
      filter(str_detect(basename(file), paste0("^", which, "_"))) %>%
      filter(str_detect(basename(file), paste0("_(?:", focal.pop, ")_")))
  } else {
    # inferred kinships depend on method -> read all
    search.dirs <- names(dir.paths.map)
    file.tbl <- tibble(dir = search.dirs) %>%
      mutate(
        kin.dir = file.path(dir, "ibd_kinship_mats"),
        file = map(kin.dir, ~ list.files(
          .x, full.names = TRUE,
          pattern = "\\.rds$"
        ))
      ) %>%
      unnest(file) %>%
      filter(!is.na(file)) %>%
      filter(!str_detect(file, "ibdendshot")) %>%
      filter(str_detect(basename(file), paste0("^", which, "_"))) %>%
      filter(str_detect(basename(file), paste0("_(?:", focal.pop, ")_")))
  }

  
  # parallel read + annotate per file
  plan(multisession, workers = cores)
  
  kin.res <- with_progress({
    p <- progressor(along = seq_len(nrow(file.tbl)))
    
    future_map_dfr(
      seq_len(nrow(file.tbl)),
      function(i) {
        p()
        
        dir   <- file.tbl$dir[i]
        f     <- file.tbl$file[i]
        fname <- basename(f)
        
        # parse: {true|inf}_{POP}_{REP}.{FILT}.rds
        m <- str_match(fname, "^(true|inf)_([^_]+)_(\\d+)\\.([^.]+)\\.rds$")
        pop  <- m[, 3]
        rep  <- as.integer(m[, 4])
        filt <- m[, 5]
        
        mat <- readRDS(f)
        
        # method metadata from your mapping
        # for true.dir that may not be in the map, guard safely:
        software.code <- unname(dir.paths.map[dir])
        if (is.na(software.code) || length(software.code) == 0) software.code <- "true"
        
        kin.tbl <- kinship.mat2tbl(mat, rep, pop)
        
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
  return(kin.res)
}
