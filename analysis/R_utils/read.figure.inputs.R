# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# read.figure.inputs.R
# ______________________________________________________________________________


library(tidyverse)
library(furrr)
library(future)

source(if (file.exists("analysis/R_utils/figure.core.R")) {
  "analysis/R_utils/figure.core.R"
} else "R_utils/figure.core.R", local = TRUE)

required.files <- function(dir, pattern) {
  files <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (!length(files)) stop("Missing required inputs in ", dir,
                           " (pattern: ", pattern, ")")
  files
}

figure.map <- function(items, fn, cores) {
  if (length(cores) != 1L || is.na(cores) || cores < 1L) cores <- 1L
  if (cores == 1L) return(purrr::map_dfr(items, fn))
  previous <- future::plan()
  on.exit(future::plan(previous), add = TRUE)
  future::plan(multisession, workers = cores)
  furrr::future_map_dfr(items, fn, .options = furrr_options(seed = TRUE))
}

read.population.files <- function(pop.info.dir, chromosome, king = FALSE) {
  stem <- paste0("^HapMapII_GRCh38_(\\d+)_chr", chromosome,
                 "_([^_]+)")
  directory <- if (king) file.path(pop.info.dir, "plink2_king_out") else {
    pop.info.dir
  }
  suffix <- if (king) "\\.kin0$" else "_unrel_samples\\.txt$"
  files <- required.files(directory, paste0(stem, suffix))
  purrr::map_dfr(files, function(file) {
    parts <- stringr::str_match(basename(file), paste0(stem, suffix))
    data <- if (king) read.plink2.king(file) else {
      tibble(ids = list(readr::read_tsv(
        file, col_names = FALSE, col_types = cols(.default = col_character()),
        show_col_types = FALSE
      )[[1]]))
    }
    data %>% mutate(rep = as.integer(parts[, 2]), pop = parts[, 3])
  })
}

read.figure.calls <- function(dir.paths.map, unrels.tbl, cores = 1L,
                              focal.pop = NULL, rep.idx = NULL) {
  populations <- unrels.tbl %>% select(pop, rep)
  if (!is.null(focal.pop)) populations <- populations %>%
    filter(pop %in% focal.pop)
  if (!is.null(rep.idx)) populations <- populations %>% filter(rep %in% rep.idx)
  if (!nrow(populations)) stop("Missing population/replicate metadata")
  datasets <- purrr::map_dfr(names(dir.paths.map), function(directory) {
    metadata <- profile.metadata(directory, dir.paths.map[[directory]], param.map)
    files <- required.files(directory,
      "_\\d+(\\.(hot|ibdends|hotibdends|ibdendshot))?\\.bed\\.gz$")
    index <- segment.metadata(files, metadata) %>%
      semi_join(populations, by = c("pop", "rep"))
    expected <- tidyr::crossing(populations,
      filt = if (metadata$param.code == "true") "none" else filt.levels)
    missing <- anti_join(expected, index, by = c("pop", "rep", "filt"))
    if (nrow(missing)) stop("Missing required segments in ", directory, ": ",
      paste(paste(missing$pop, missing$rep, missing$filt, sep = "/"),
            collapse = ", "))
    index
  })
  result <- figure.map(seq_len(nrow(datasets)), function(i) {
    metadata <- datasets[i, ]
    ids <- unrels.tbl$ids[unrels.tbl$rep == metadata$rep &
                           unrels.tbl$pop == metadata$pop]
    if (length(ids) != 1L) stop("Missing or duplicate unrelated IDs: ",
                                metadata$pop, " replicate ", metadata$rep)
    data <- readr::read_tsv(metadata$file, show_col_types = FALSE,
      col_types = cols(.default = col_double(), chr = col_character(),
                       sam1 = col_character(), sam2 = col_character()))
    bind_cols(data, metadata[rep(1L, nrow(data)), setdiff(names(metadata), "file")]) %>%
      mutate(is.unrel = sam1 %in% ids[[1]] & sam2 %in% ids[[1]])
  }, cores)
  attr(result, "datasets") <- datasets
  result
}

read.figure.kinship <- function(dir.paths.map, cores = 1L,
                                which = "inf", focal.pop = NULL,
                                rep.idx = NULL, true.dir = NULL) {
  if (which == "true") {
    directory <- if (is.null(true.dir)) names(dir.paths.map)[1] else true.dir
    dir.paths.map <- setNames("true", directory)
  }
  datasets <- purrr::map_dfr(names(dir.paths.map), function(directory) {
    metadata <- profile.metadata(directory, dir.paths.map[[directory]], param.map)
    files <- required.files(file.path(directory, "ibd_kinship_mats"),
      paste0("^", which, "_[^_]+_\\d+\\.[^.]+\\.rds$"))
    parts <- str_match(basename(files),
      "^(true|inf)_([^_]+)_(\\d+)\\.([^.]+)\\.rds$")
    tibble(file = files, directory = directory, pop = parts[, 3],
           rep = as.integer(parts[, 4]), filt = parts[, 5]) %>%
      bind_cols(metadata[rep(1L, length(files)), ])
  })
  if (!is.null(focal.pop)) datasets <- datasets %>% filter(pop %in% focal.pop)
  if (!is.null(rep.idx)) datasets <- datasets %>% filter(rep %in% rep.idx)
  if (!nrow(datasets)) stop("Missing required kinship population/replicate")
  axes <- datasets %>% distinct(pop, rep)
  expected <- tidyr::crossing(axes, directory = names(dir.paths.map),
                             filt = if (which == "true") "none" else filt.levels)
  missing <- anti_join(expected, datasets,
                       by = c("directory", "pop", "rep", "filt"))
  if (nrow(missing)) stop("Missing required kinship input: ",
    file.path(missing$directory[1], "ibd_kinship_mats",
      paste0(which, "_", missing$pop[1], "_", missing$rep[1], ".",
             missing$filt[1], ".rds")))
  if (which == "true") datasets <- datasets %>% filter(filt == "none")
  figure.map(seq_len(nrow(datasets)), function(i) {
    metadata <- datasets[i, ]
    data <- kinship.mat2tbl(readRDS(metadata$file), metadata$rep, metadata$pop)
    bind_cols(data, metadata[rep(1L, nrow(data)),
        c("filt", "software.code", "param.code", "software.with.params", "imap")])
  }, cores)
}
